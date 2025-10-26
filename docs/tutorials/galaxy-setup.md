# Galaxy Setup Tutorial

Learn how to set up a local Galaxy instance with COBRAxy tools for web-based metabolic analysis.

## Overview

This tutorial provides guidance and references to official Galaxy documentation for:

- Installing Galaxy locally (using official Galaxy guides)
- Adding COBRAxy tools to your Galaxy instance
- Running COBRAxy analyses through the web interface

## Step 1: Install Galaxy Locally

For installing Galaxy on your local machine, follow the official documentation:

### Official Installation Guides

- **[Galaxy Installation Guide](https://docs.galaxyproject.org/en/master/admin/)**

- **[Galaxy Quick Start](https://docs.galaxyproject.org/en/master/admin/production.html)**

- **[Galaxy Training: Admin Track](https://training.galaxyproject.org/training-material/topics/admin/)**

### Quick Summary

1. Clone Galaxy repository: `git clone https://github.com/galaxyproject/galaxy.git`
2. Run setup script: `sh run.sh`  
3. Access at: `http://localhost:8080`

**Note**: Refer to official documentation for detailed configuration, security settings, and production deployment.

## Step 2: Add COBRAxy Tools to Galaxy

For adding custom tools to Galaxy, refer to the official documentation:

### Official Tool Installation Guides

- **[Adding Tools to Galaxy](https://docs.galaxyproject.org/en/master/admin/tool_panel.html)**

- **[Tool Development Tutorial](https://training.galaxyproject.org/training-material/topics/dev/)**

- **[Galaxy Tool Installation](https://planemo.readthedocs.io/en/latest/)**

### COBRAxy-Specific Setup

1. **Link COBRAxy to Galaxy tools** directory:
   ```bash
   cd /path/to/galaxy
   ln -s /path/to/COBRAxy/src tools/cobraxy
   ```

2. **Add tools to Galaxy configuration**:
   Edit `config/tool_conf.xml` and add a COBRAxy section:
   ```xml
   <section id="cobraxy" name="COBRAxy">
     <tool file="cobraxy/importMetabolicModel.xml" />
     <tool file="cobraxy/exportMetabolicModel.xml" />
     <tool file="cobraxy/ras_generator.xml" />
     <tool file="cobraxy/rps_generator.xml" />
     <tool file="cobraxy/marea.xml" />
     <tool file="cobraxy/ras_to_bounds.xml" />
     <tool file="cobraxy/flux_simulation.xml" />
     <tool file="cobraxy/flux_to_map.xml" />
     <tool file="cobraxy/marea_cluster.xml" />
   </section>
   ```

3. **Restart Galaxy** to load the new tools.

**Note**: Consult the official Galaxy documentation for detailed instructions on tool installation, dependency management, and troubleshooting.

## Step 3: Using COBRAxy in Galaxy

### Verify Installation

After following the official Galaxy setup and tool installation procedures:

1. Access your Galaxy instance (typically `http://localhost:8080`)
2. Check that COBRAxy tools appear in the tool panel
3. Verify Python dependencies are available in Galaxy's environment

### Basic Usage

1. **Upload data** using Galaxy's data upload interface
2. **Select COBRAxy tools** from the tool panel
3. **Configure parameters** through the web interface
4. **Execute analyses** and monitor job progress
5. **Download results** from Galaxy's history panel

## Creating COBRAxy Workflows

### Workflow Development Resources

For creating workflows with COBRAxy tools in Galaxy:

- **[Galaxy Workflow Tutorial](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/workflow-editor/tutorial.html)**
  - Creating, editing, and sharing workflows
  - Workflow best practices

- **[Workflow Management](https://docs.galaxyproject.org/en/master/user/galaxy_workflow.html)**
  - Official workflow documentation
  - Advanced workflow features

### Example COBRAxy Workflow

A typical COBRAxy workflow might include:

1. **RAS Generator** → Generate activity scores from gene expression
2. **MAREA** → Perform statistical analysis and create pathway maps
3. **RAS to Bounds** → Apply constraints (optional, for flux analysis)
4. **Flux Simulation** → Sample metabolic fluxes (optional)
5. **Flux to Map** → Create final visualizations (optional)

## Additional Resources

### Galaxy Administration Resources

- **[Galaxy Admin Documentation](https://docs.galaxyproject.org/en/master/admin/)**
  - Complete administrator guide
  - Configuration, security, and maintenance

- **[Galaxy Training Materials](https://training.galaxyproject.org/)**
  - Hands-on tutorials for administrators and users
  - Best practices and troubleshooting

- **[Galaxy Community Hub](https://galaxyproject.org/)**
  - Community support and resources
  - Tool repositories and shared workflows

### COBRAxy-Specific Resources

- **Dependencies**: Ensure `cobra`, `pandas`, `numpy`, `scipy` are installed in Galaxy's Python environment
- **Tool Files**: All COBRAxy XML and Python files should be accessible to Galaxy
- **Configuration**: Follow Galaxy's tool installation procedures for proper integration

## Troubleshooting

For troubleshooting Galaxy installations and tool integration issues:

### Official Troubleshooting Resources

- **[Galaxy FAQ](https://docs.galaxyproject.org/en/master/admin/faq.html)**
  - Common installation and configuration issues
  - Performance optimization tips

- **[Galaxy Help Forum](https://help.galaxyproject.org/)**
  - Community-driven support
  - Search existing solutions or ask new questions

- **[Galaxy GitHub Issues](https://github.com/galaxyproject/galaxy/issues)**
  - Report bugs and technical issues
  - Feature requests and discussions

### COBRAxy-Specific Issues

For issues specific to COBRAxy tools in Galaxy:

- **Tool not appearing**: Check tool_conf.xml configuration and restart Galaxy
- **Execution failures**: Verify Python dependencies and file permissions  
- **Parameter errors**: Ensure input data formats match tool requirements

Refer to the [COBRAxy Tools Documentation](/tools/) for detailed parameter information and data format requirements.

## Summary

This tutorial provides guidance for setting up Galaxy with COBRAxy tools by referencing official Galaxy documentation. For detailed installation procedures, always consult the official Galaxy administrator guides, as they are regularly updated with the latest best practices and troubleshooting information.

The combination of Galaxy's web interface with COBRAxy's metabolic analysis capabilities provides a powerful platform for researchers who prefer graphical interfaces over command-line tools.