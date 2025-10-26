from setuptools import setup, find_packages
import os

# Get the path to README.md in the parent directory
readme_path = os.path.join(os.path.dirname(__file__), '..', 'README.md')

setup(
    name='cobraxy',
    version='0.1.0',
    description='A collection of tools for metabolic flux analysis in Galaxy.',
    long_description=open(readme_path, encoding="utf-8").read(),
    long_description_content_type='text/markdown',
    author='',  
    author_email='',
    url='https://github.com/CompBtBs/COBRAxy.git',
    license='',
    packages=find_packages(include=["utils", "utils.*"]),  
    py_modules=[
        'importMetabolicModel',
        'exportMetabolicModel',
        'ras_generator',
        'rps_generator',
        'marea_cluster',
        'marea',
        'ras_to_bounds',
        'flux_simulation',
        'flux_to_map'
    ],
    include_package_data=True, 
    install_requires=[
        'cairosvg>=2.7.0',
        'cobra>=0.29.0',
        'joblib>=1.3.0',
        'lxml>=5.0.0',
        'matplotlib>=3.7.0',
        'numpy>=1.24.0',
        'pandas>=2.0.0',
        'pyvips>=2.2.0',
        'scikit-learn>=1.3.0',
        'scipy>=1.11.0',
        'seaborn>=0.13.0',
        'svglib>=1.5.0',
        'anndata>=0.8.0',
        'pydeseq2>=0.4.0'
    ],
    entry_points={
        'console_scripts': [
            'importMetabolicModel=importMetabolicModel:main',
            'exportMetabolicModel=exportMetabolicModel:main',
            'ras_generator=ras_generator:main',
            'rps_generator=rps_generator:main',
            'marea_cluster=marea_cluster:main',
            'marea=marea:main',
            'ras_to_bounds=ras_to_bounds:main',
            'flux_simulation=flux_simulation:main',
            'flux_to_map=flux_to_map:main'
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8,<3.14',
)
