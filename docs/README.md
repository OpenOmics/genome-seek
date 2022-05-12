# Build documentation  

> **Please Note:** When a commit is pushed to the `docs/` directory, it triggers a [github actions workflow](https://github.com/OpenOmics/genome-seek/actions) to build the static-site and push it to the gh-pages branch.

### Installation
```bash
# Clone the Repository
git clone https://github.com/OpenOmics/genome-seek.git
# Create a virtual environment
python3 -m venv .venv
# Activate the virtual environment
. .venv/bin/activate
# Update pip
pip install --upgrade pip
# Download Dependencies
pip install -r docs/requirements.txt
```

### Preview while editing  
MkDocs includes a previewing server, so you can view your updates live and as you write your documentation. The server will automatically rebuild the site upon editing and saving a file.  
```bash
# Activate the virtual environment
. .venv/bin/activate
# Start serving your documentation
mkdocs serve
```

### Build static site  
Once you are content with your changes, you can build the static site:  
```bash
mkdocs build
```
