from setuptools import setup, find_packages

setup(
    name='cobraxy',
    version='0.1.0',
    description='A collection of tools for metabolic flux analysis in Galaxy.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='',  
    author_email='',
    url='https://github.com/CompBtBs/COBRAxy.git',
    license='',
    packages=find_packages(include=["utils", "utils.*"]),  
    py_modules=[
        'ras_generator',
        'rps_generator',
        'marea_cluster',
        'marea',
        'metabolic_model_setting',
        'ras_to_bounds',
        'flux_simulation',
        'flux_to_map'
    ],
    include_package_data=True, 
    install_requires=[
        'cairosvg==2.7.1',
        'cobra==0.29.0',
        'joblib==1.4.2',
        'lxml==5.2.2',
        'matplotlib==3.7.3',
        'numpy==1.24.4',
        'pandas==2.0.3',
        'pyvips==2.2.3',
        'scikit-learn==1.3.2',
        'scipy==1.11',
        'seaborn==0.13.0',
        'svglib==1.5.1',
        'anndata==0.8.0',
        'pydeseq2==0.5.1'
    ],
    entry_points={
        'console_scripts': [
            'metabolic_model_setting=metabolic_model_setting:main',
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
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8.20,<3.12',
)
