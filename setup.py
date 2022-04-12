import setuptools

setuptools.setup(
    name='gisaid_script',
    version='0.1',
    description='A package for formatting metadata for SARS-CoV-2 metadata',
    long_description=(
        ''''''
    ),
    long_description_content_type='text/markdown',
    url='https://github.com/mackellardrew/gisaid_script',
    author='Drew MacKellar',
    author_email='drew.mackellar@doh.wa.gov',
    license='GNU GPL',
    # packages=['gisaid_script'],
    packages=setuptools.find_packages(exclude=['docs', 'tests*']),
    install_requires=[
        'pandas', 'gsutil', 'biopython', 'IPython', 'openpyxl',
    ],
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    extras_require = {
        'dev': [
            'pytest>=3.7',
        ],
    },
    # scripts=['bin/gisaid_script'],
    entry_points={
        'console_scripts': [
            'gisaid_script=gisaid_script.command_line:cmd_line_main',
        ],
    },
    zip_safe=False
)