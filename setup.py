import setuptools

setuptools.setup(
    name='gisaid_script',
    version='0.1',
    description='A package for formatting metadata for SARS-CoV-2 metadata',
    url='https://github.com/mackellardrew/gisaid_script',
    author='Drew MacKellar',
    author_email='drew.mackellar@doh.wa.gov',
    license='GNU GPL',
    packages=['gisaid_script'],
    install_requires=[
        'pandas', 'gsutil', 'biopython', 'IPython'
    ],
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    scripts=['bin/gisaid_script'],
    zip_safe=False
)