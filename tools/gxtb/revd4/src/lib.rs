// SPDX-License-Identifier: LGPL-3.0-or-later
//! Periodic image sums for the g-xTB 2.0.1 revD4 model.
use axtb_base::Molecule;
use axtb_gxtb::model::Model;
use axtb_gxtb::options::Conventions;
use axtb_gxtb::revd4::{self, RevD4};
use std::ffi::{c_char, c_int};
use std::{ptr, slice};

#[cfg(test)]
mod tests;

type V = [f64; 3];
fn dot(a: V, b: V) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

struct Pair {
    i: usize,
    j: usize,
    d: V,
    r: f64,
    dcn: f64,
    dk_i: f64,
    dk_j: f64,
    radial: f64,
}

struct Derivatives {
    grad: Vec<V>,
    strain: [V; 3],
}
impl Derivatives {
    fn new(n: usize) -> Self {
        Self {
            grad: vec![[0.; 3]; n],
            strain: [[0.; 3]; 3],
        }
    }
    fn add_image(&mut self, atom: usize, g: V, d: V) {
        for a in 0..3 {
            self.grad[atom][a] += g[a];
            for b in 0..3 {
                self.strain[b][a] += g[a] * d[b];
            }
        }
    }
    fn add_pair(&mut self, p: &Pair, dedr: f64) {
        let g = p.d.map(|v| v * dedr / p.r);
        self.add_image(p.j, g, p.d);
        self.add_image(p.i, g.map(|v| -v), [0.; 3]);
    }
}

pub struct Geometry {
    z: Vec<usize>,
    pairs: Vec<Pair>,
    weights: Vec<f64>,
    atm: Vec<f64>,
    atm_derivatives: Derivatives,
}

impl Geometry {
    pub fn new(z: &[usize], xyz: &[V], cell: [V; 3], pbc: [bool; 3]) -> Result<Self, String> {
        let n = z.len();
        if n == 0
            || !z.iter().all(|x| (1..=103).contains(x))
            || xyz.len() != n
            || !xyz
                .iter()
                .flatten()
                .chain(cell.iter().flatten())
                .all(|v| v.is_finite())
        {
            return Err("invalid revD4 geometry or element outside Z=1..103".into());
        }
        let model = Model::v2_0_1(Conventions::Gxtb);
        let mut reciprocal = [[0.; 3]; 3];
        if pbc.contains(&true) {
            let volume = dot(cell[0], cross(cell[1], cell[2]));
            if !volume.is_finite() || volume.abs() < 1e-12 {
                return Err("singular periodic cell".into());
            }
            for a in 0..3 {
                reciprocal[a] = cross(cell[(a + 1) % 3], cell[(a + 2) % 3]).map(|v| v / volume);
            }
        }
        let mut cn = vec![0.; n];
        let mut pairs = Vec::new();
        let mut near: Vec<Vec<(usize, V)>> = (0..n).map(|i| vec![(i, [0.; 3])]).collect();
        for i in 0..n {
            for j in i..n {
                let delta: V = std::array::from_fn(|a| xyz[j][a] - xyz[i][a]);
                let mut lo = [0i32; 3];
                let mut hi = [0i32; 3];
                for a in 0..3 {
                    if pbc[a] {
                        let center = -dot(delta, reciprocal[a]);
                        let reach = revd4::CUTOFF_2B * dot(reciprocal[a], reciprocal[a]).sqrt();
                        lo[a] = (center - reach - 1e-12).ceil() as i32;
                        hi[a] = (center + reach + 1e-12).floor() as i32;
                    }
                }
                for u in lo[0]..=hi[0] {
                    for v in lo[1]..=hi[1] {
                        for w in lo[2]..=hi[2] {
                            if i == j && [u, v, w] <= [0, 0, 0] {
                                continue;
                            }
                            let d: V = std::array::from_fn(|a| {
                                delta[a]
                                    + f64::from(u) * cell[0][a]
                                    + f64::from(v) * cell[1][a]
                                    + f64::from(w) * cell[2][a]
                            });
                            let r = dot(d, d).sqrt();
                            if r > revd4::CUTOFF_2B {
                                continue;
                            }
                            if r < 1e-6 {
                                return Err("coincident atoms or periodic images".into());
                            }
                            let mut dcn = 0.;
                            if r <= revd4::CUTOFF_CN {
                                let rc = model.covalent_rad_d3(z[i]) + model.covalent_rad_d3(z[j]);
                                let count = axtb_d4::cn_count(rc, r);
                                cn[i] += count;
                                cn[j] += count;
                                dcn = axtb_d4::cn_dcount(rc, r);
                            }
                            if r <= revd4::CUTOFF_3B {
                                near[i].push((j, d));
                                near[j].push((i, d.map(|x| -x)));
                            }
                            pairs.push(Pair {
                                i,
                                j,
                                d,
                                r,
                                dcn,
                                dk_i: 0.,
                                dk_j: 0.,
                                radial: 0.,
                            });
                        }
                    }
                }
            }
        }
        // The molecular API omits diagonal atom pairs. A second copy with
        // the same supplied CN exposes the coefficient between identical
        // images, including its two independent CN derivatives.
        let twin_z = [z, z].concat();
        let dummy: Vec<V> = (0..2 * n).map(|i| [2. * i as f64, 0., 0.]).collect();
        let twins = Molecule::new(&twin_z, &dummy, 0., None).map_err(|e| e.to_string())?;
        let data = RevD4::with_cn(&model, &twins, cn.repeat(2)).map_err(|e| e.to_string())?;
        let mut weights = vec![0.; n * n];
        for p in &mut pairs {
            let (i, j) = (p.i, p.j);
            let ij = i * (2 * n) + n + j;
            let ji = j * (2 * n) + n + i;
            let rho = (3. * data.r4r2[i] * data.r4r2[j]).sqrt();
            let (reff, dr) = revd4::reff(p.r, rho);
            let inv = 1. / (p.r + reff);
            let damp = revd4::S6 * inv.powi(6);
            let weight = data.kmat[ij] * damp;
            weights[i * n + j] += weight;
            weights[j * n + i] += weight;
            p.dk_i = -damp * data.dkmat[ij];
            p.dk_j = -damp * data.dkmat[ji];
            p.radial = 6. * weight * inv * (1. + dr);
        }
        let mut state = Self {
            z: z.to_vec(),
            pairs,
            weights,
            atm: vec![0.; n],
            atm_derivatives: Derivatives::new(n),
        };
        let mut dcn = vec![0.; n];
        for (i, images) in near.iter().enumerate() {
            if images.len() < 3 {
                continue;
            }
            let ids: Vec<_> = images.iter().map(|(j, _)| z[*j]).collect();
            let pos: Vec<_> = images.iter().map(|(_, d)| *d).collect();
            let env = Molecule::new(&ids, &pos, 0., None).map_err(|e| e.to_string())?;
            let local_cn = images.iter().map(|(j, _)| cn[*j]).collect();
            let kernel = RevD4::with_cn(&model, &env, local_cn).map_err(|e| e.to_string())?;
            let mut grad = vec![[0.; 3]; images.len()];
            let (e, dc) = kernel.three_body_of(&model, &env, Some(&mut grad), Some(0));
            // Each physical triple is anchored at its three vertices.
            state.atm[i] = e / 3.;
            for (k, &(j, d)) in images.iter().enumerate() {
                dcn[j] += dc[k] / 3.;
                state
                    .atm_derivatives
                    .add_image(j, grad[k].map(|v| v / 3.), d);
            }
        }
        for p in &state.pairs {
            state
                .atm_derivatives
                .add_pair(p, p.dcn * (dcn[p.i] + dcn[p.j]));
        }
        Ok(state)
    }

    pub fn electronic(&self, q: &[f64], energy: &mut [f64], potential: &mut [f64]) {
        let n = self.z.len();
        let f: Vec<_> = self
            .z
            .iter()
            .zip(q)
            .map(|(&z, &q)| revd4::charge_function(z, q))
            .collect();
        for (i, row) in self.weights.chunks(n).enumerate() {
            let sum: f64 = row.iter().zip(&f).map(|(w, (v, _))| w * v).sum();
            energy[i] = -0.5 * f[i].0 * sum;
            potential[i] = -f[i].1 * sum;
        }
    }

    fn electronic_derivatives(&self, q: &[f64]) -> Derivatives {
        let mut out = Derivatives::new(self.z.len());
        let f: Vec<_> = self
            .z
            .iter()
            .zip(q)
            .map(|(&z, &q)| revd4::charge_function(z, q).0)
            .collect();
        let mut dcn = vec![0.; self.z.len()];
        for p in &self.pairs {
            let factor = f[p.i] * f[p.j];
            dcn[p.i] += factor * p.dk_i;
            dcn[p.j] += factor * p.dk_j;
            out.add_pair(p, factor * p.radial);
        }
        for p in &self.pairs {
            out.add_pair(p, p.dcn * (dcn[p.i] + dcn[p.j]));
        }
        out
    }
}

// Fortran owns one opaque Box per container cache. Buffers use atomic
// units and contiguous Fortran layout xyz(3,n), lattice(3,3), sigma(3,3).
#[no_mangle]
pub unsafe extern "C" fn cp2k_revd4_new(
    n: c_int,
    z: *const c_int,
    xyz: *const V,
    cell: *const [V; 3],
    pbc: *const [c_int; 3],
    msg: *mut c_char,
    len: usize,
) -> *mut Geometry {
    let result = std::panic::catch_unwind(|| {
        if n <= 0 {
            return Err("empty revD4 input".to_string());
        }
        let numbers: Vec<_> = slice::from_raw_parts(z, n as usize)
            .iter()
            .map(|&v| v as usize)
            .collect();
        Geometry::new(
            &numbers,
            slice::from_raw_parts(xyz, n as usize),
            *cell,
            (*pbc).map(|x| x != 0),
        )
    })
    .unwrap_or_else(|_| Err("revD4 geometry preparation failed".into()));
    match result {
        Ok(state) => Box::into_raw(Box::new(state)),
        Err(err) => {
            if len > 0 {
                let count = err.len().min(len - 1);
                ptr::copy_nonoverlapping(err.as_ptr(), msg.cast(), count);
                *msg.add(count) = 0;
            }
            ptr::null_mut()
        }
    }
}
#[no_mangle]
pub unsafe extern "C" fn cp2k_revd4_free(state: *mut Geometry) {
    if !state.is_null() {
        drop(Box::from_raw(state));
    }
}
#[no_mangle]
pub unsafe extern "C" fn cp2k_revd4_scf(state: &Geometry, q: *const f64, e: *mut f64, v: *mut f64) {
    let n = state.z.len();
    state.electronic(
        slice::from_raw_parts(q, n),
        slice::from_raw_parts_mut(e, n),
        slice::from_raw_parts_mut(v, n),
    );
}
#[no_mangle]
pub unsafe extern "C" fn cp2k_revd4_derivatives(
    state: &Geometry,
    q: *const f64,
    g: *mut V,
    s: *mut [V; 3],
) {
    let n = state.z.len();
    let d = state.electronic_derivatives(slice::from_raw_parts(q, n));
    ptr::copy_nonoverlapping(d.grad.as_ptr(), g, n);
    *s = d.strain;
}
#[no_mangle]
pub unsafe extern "C" fn cp2k_revd4_atm(state: &Geometry, e: *mut f64, g: *mut V, s: *mut [V; 3]) {
    let n = state.z.len();
    ptr::copy_nonoverlapping(state.atm.as_ptr(), e, n);
    ptr::copy_nonoverlapping(state.atm_derivatives.grad.as_ptr(), g, n);
    *s = state.atm_derivatives.strain;
}
