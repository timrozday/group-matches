from setuptools import setup

setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='group_matches',
    url='https://github.com/timrozday/group_matches',
    author='Tim Rozday',
    author_email='timrozday@ebi.ac.uk',
    # Needed to actually package something
    packages=['group_matches'],
    # Needed for dependencies
    install_requires=[],
    version='0.1',
    # The license can be anything you like
    license='Do what you like with it (just nothing evil)',
    description='A few functions for grouping matches from text_query.',
    # We will also need a readme eventually (there will be a warning)
    long_description='A few functions for grouping matches from text_query.',
    # long_description=open('README.txt').read(),
)
