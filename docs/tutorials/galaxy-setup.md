# Galaxy Setup Tutorial

Learn how to set up a local Galaxy instance with COBRAxy tools for web-based metabolic analysis.

## Overview

Galaxy provides a web interface for running COBRAxy tools without command line knowledge. This tutorial covers:

- Installing Galaxy locally
- Adding COBRAxy tools to Galaxy
- Running your first analysis through the web interface
- Creating reusable workflows

**Time required**: ~30 minutes  
**Difficulty**: Beginner  
**Prerequisites**: Basic computer administration skills

## Step 1: Install Galaxy

### Download and Setup Galaxy

```bash
# Clone Galaxy repository (stable release)
git clone -b release_23.1 https://github.com/galaxyproject/galaxy.git
cd galaxy

# Copy sample configuration
cp config/galaxy.yml.sample config/galaxy.yml
```

### Configure Galaxy

Edit `config/galaxy.yml` to customize your installation:

```yaml
galaxy:
  # Basic configuration
  admin_users: admin@galaxy.local
  database_connection: sqlite:///database/universe.sqlite?isolation_level=IMMEDIATE
  file_path: database/files
  tool_dependency_dir: database/dependencies
  
  # Security settings  
  id_secret: your-secret-key-here
  use_remote_user: false
  
  # Job configuration
  job_config_file: config/job_conf.yml
```

### Start Galaxy

```bash
# Start Galaxy server
sh run.sh

# Galaxy will be available at: http://localhost:8080
# First startup takes 10-15 minutes to initialize database
```

**Note**: Keep the terminal open - Galaxy runs in the foreground.

## Step 2: Install COBRAxy Tools

### Method 1: Copy Tool Files

```bash
# Create COBRAxy tool directory
mkdir -p tools/cobraxy

# Copy COBRAxy XML files (adjust path to your COBRAxy installation)
cp /path/to/COBRAxy/*.xml tools/cobraxy/
cp -r /path/to/COBRAxy/utils tools/cobraxy/
cp -r /path/to/COBRAxy/local tools/cobraxy/
```

### Method 2: Symlink (Development)

```bash
# Create symlinks for development
mkdir -p tools/cobraxy
ln -s /path/to/COBRAxy/*.xml tools/cobraxy/
ln -s /path/to/COBRAxy/*.py tools/cobraxy/
ln -s /path/to/COBRAxy/utils tools/cobraxy/utils
ln -s /path/to/COBRAxy/local tools/cobraxy/local
```

### Register Tools in Galaxy

Edit `config/tool_conf.xml` to add COBRAxy tools:

```xml
<?xml version='1.0' encoding='utf-8'?>
<toolbox monitor="true">
  <!-- Other tool sections... -->
  
  <section id="cobraxy" name="COBRAxy - Metabolic Analysis">
    <tool file="cobraxy/ras_generator.xml" />
    <tool file="cobraxy/rps_generator.xml" />
    <tool file="cobraxy/marea.xml" />
    <tool file="cobraxy/ras_to_bounds.xml" />
    <tool file="cobraxy/flux_simulation.xml" />
    <tool file="cobraxy/flux_to_map.xml" />
    <tool file="cobraxy/metabolicModel2Tabular.xml" />
    <tool file="cobraxy/tabular2MetabolicModel.xml" />
    <tool file="cobraxy/marea_cluster.xml" />
  </section>
</toolbox>
```

### Restart Galaxy

```bash
# Stop Galaxy (Ctrl+C in terminal)
# Restart
sh run.sh
```

## Step 3: Verify Installation

### Access Galaxy Interface

1. Open http://localhost:8080 in your browser
2. Create an admin account (use email from `admin_users` config)
3. Look for "COBRAxy - Metabolic Analysis" in the left tool panel

### Test Tool Installation

1. Click on **RAS Generator** in the tool panel
2. You should see the tool interface with parameter options
3. If you see error messages, check:
   - File paths in tool XMLs are correct
   - Python dependencies are installed
   - COBRAxy files have proper permissions

## Step 4: Upload Sample Data

### Prepare Test Dataset

Create a sample gene expression file:

```tsv
Gene_ID	Control_1	Control_2	Treatment_1	Treatment_2
HGNC:5	10.5	11.2	15.7	14.3
HGNC:10	3.2	4.1	8.8	7.9
HGNC:15	7.9	8.2	4.4	5.1
HGNC:25	12.1	13.5	18.2	17.8
```

### Upload to Galaxy

1. Click **Upload Data** icon (📁) in Galaxy
2. Choose **Choose local files**
3. Select your test file
4. Set **Type** to `tabular`
5. Click **Start** to upload

## Step 5: Run Your First Analysis

### Generate RAS Scores

1. Select **RAS Generator** from tools panel
2. Configure parameters:
   - **Input dataset**: Your uploaded gene expression file
   - **Rule selector**: ENGRO2
   - **Output name**: `ras_scores`
3. Click **Execute**

### Monitor Job Progress

- Jobs appear in the **History** panel (right side)
- Wait for job to turn green (completed)
- Red indicates error - check job details

### Generate Pathway Maps

1. Select **MAREA** from tools panel
2. Configure parameters:
   - **Use RAS**: Yes
   - **Input data**: Select your RAS output
   - **Map choice**: ENGRO2
   - **Gene set analysis**: Yes
3. Click **Execute**

## Step 6: Create Reusable Workflows

### Extract Workflow

1. Go to **Workflow** menu → **Extract Workflow**
2. Select history containing your analysis
3. Name workflow: "COBRAxy Basic Analysis"
4. Click **Create Workflow**

### Edit Workflow

1. Go to **Workflow** menu → **Edit**
2. Select your workflow
3. Add tool descriptions and annotations
4. Connect tools with clear parameter flow
5. Save workflow

### Share Workflow

1. Go to **Workflow** menu → **Share or Publish**
2. Make workflow public or share with collaborators
3. Export workflow file for distribution

## Advanced Configuration

### Install Dependencies

Install required packages in Galaxy's environment:

```bash
# Install in Galaxy's Python environment
source .venv/bin/activate  # or your Galaxy venv

# Install COBRAxy dependencies
pip install cobra pandas numpy scipy

# Install optional solver
pip install swiglpk  # GLPK solver
```

### Configure Job Runners

Edit `config/job_conf.yml` for better job management:

```yaml
runners:
  local:
    load: galaxy.jobs.runners.local:LocalJobRunner
    workers: 4

execution:
  default: local
  environments:
    local:
      runner: local
      tmp_dir: true
```

### Set Up Tool Dependencies

```bash
# Install via conda (recommended)
conda install -c conda-forge cobra glpk

# Or use Galaxy's dependency system
mkdir -p database/dependencies
# Install tools through Galaxy admin interface
```

## Troubleshooting

### Common Issues

**COBRAxy tools not appearing**
- Check `tool_conf.xml` syntax
- Verify file paths in XML files
- Restart Galaxy after changes

**Tool execution failures**
- Check Galaxy logs: `tail -f main.log`
- Verify Python dependencies installed
- Check file permissions on COBRAxy files

**Slow execution**
- Increase worker count in `job_conf.yml`
- Install optional dependencies (GLPK, etc.)
- Check system resources

### Debug Mode

Run Galaxy in debug mode for detailed error messages:

```bash
# Enable debug mode
export GALAXY_CONFIG_DEBUG=true
sh run.sh
```

### Log Files

Check these log files for troubleshooting:
- `main.log` - Galaxy server logs
- `handler.log` - Job execution logs
- `uwsgi.log` - Web server logs

## Next Steps

Now that Galaxy is set up:

1. **[Basic Workflow Tutorial](workflow.md)** - Learn complete analysis pipeline
2. **[Data Formats Guide](data-formats.md)** - Prepare your own data
3. **Create custom workflows** for your specific analysis needs
4. **Set up multi-user access** for your lab or organization

## Resources

- [Galaxy Documentation](https://docs.galaxyproject.org/)
- [Galaxy Training Materials](https://training.galaxyproject.org/)
- [Galaxy Tool Development](https://planemo.readthedocs.io/)
- [COBRAxy Galaxy Tools](https://github.com/CompBtBs/COBRAxy/tree/main/Galaxy_tools)