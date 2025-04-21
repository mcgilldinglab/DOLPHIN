# Usage

(installation)=
## Installation

**Note:** This tool has been primarily tested on Linux-based systems. While it may run on other platforms, we recommend using a Linux environment for best compatibility and performance—especially for memory-intensive preprocessing steps such as STAR or Cell Ranger alignment.

Install from Github

📥 Step 1: Clone the Repository
```bash
git clone https://github.com/mcgilldinglab/DOLPHIN.git
cd DOLPHIN
```

🛠 Step 2: Create and Activate the Conda Environment
```bash
conda env create -f environment.yaml
conda activate DOLPHIN
```

🧑‍💻 (Optional) Step 3: Developer Mode Installation
```bash
pip install -e .
```

✅ Step 4: Validate the Installation
You can check if the package is correctly installed by opening Python and running:
```bash
import DOLPHIN
```

## Example Dataset
You can download the processed dataset from [One Drive.](https://mcgill-my.sharepoint.com/:f:/g/personal/kailu_song_mail_mcgill_ca/EvZtHeW7qjJJs_RHc2-327ABeLXafa-ruvfk9Vs134crig?e=VBn7KG)
and follow the [example tutorial.](examples/run_DOLPHIN)

**Required files**

For human alignment to the exon level, please download the pre-processed gtf file here:[One Drive.](https://outlook.office.com/host/377c982d-9686-450e-9a7c-22aeaf1bc162/7211f19f-262a-42eb-a02e-289956491741)


