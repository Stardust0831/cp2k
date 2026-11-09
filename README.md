# SAI cuSOLVERMp tester

Computing resources are provided by [OPEN-SAI](https://www.open-sai.com/?lang=en-US) (Open Source Supercomputing Center of SAI).

Four V100 GPUs, four MPI ranks and four OpenMP threads per rank. Tests run inside the complete accepted daily CP2K Apptainer image.

[Latest report](report.txt) · [Detailed results](summary.md) · [Latest Action](https://github.com/Stardust0831/cp2k/actions/runs/42555942513)

Raw inputs, outputs and provenance are retained in Action artifacts for 90 days and on SAI. The source commit in the report identifies the tested image.

An isolated self-hosted runner on the relay VPS checks hourly for newly published development/16V100 images. Each image is tested once unless manually requested. SAI accepts only the fixed tester commands through a dedicated restricted SSH key.
