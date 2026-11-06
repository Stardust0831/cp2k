# Baseline full-spectrum fixture diagnostic

Job 2234785 compared the original DZVP-GTH Si8 basis with DZVP-MOLOPT-SR-GTH.
Both have 104 basis functions. Inputs, raw CP2K output, overlap-condition output,
GPU/rank records and accounting are retained. Restart wavefunctions are omitted.

The reported overlap condition number for the k-point input is for Gamma.
The real and complex MOLOPT-SR full-spectrum solves passed without warnings;
the original DZVP-GTH full spectrum showed orthogonality warnings. See the
main report for the final fixture choice and unchanged regression tolerances.
