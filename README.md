### Install
1. Clone the repo.
2. Get submodules within this repo. Read through "Cloning a Project with Submodules" in [here](https://git-scm.com/book/en/v2/Git-Tools-Submodules) if unfamiliar with submodules.
    - Initialize local submodule config - `git submodule init`
    - Fetch data from the repo and check out the appropriate commit - `git submodule update`
3. Create the conda environment - `conda env create -f environment.yml`
4. Activate conda environment - `conda activate preprocessing`

### TODOs
We originally were trying to go "top-down" by thinking of all possible functionalities needed, then finding the quickest way to implement each, then using those for specific examples we have. While not incorrect, it seems more of a challenge due to our inexperience with these problem. We'll instead take a "bottom-up" approach, where we're going to work towards specific examples which will help us realize the set of need functionalities and how to quickly implement them. See [TODO.md](./TODO.md) for more.

### Notes
Currently have have a disconnected set of repos for CMAAS TA3. Ultimately, we'll want to connect [most] of them. Should keep that in the back of our minds when developing.