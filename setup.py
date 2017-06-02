from setuptools import setup

# Get version from xyalign/_version.py

version_py = "xyalign/_version.py"
exec(open(version_py).read())

package_description = 'Tools to infer ploidy, correct for sex chromosome "\
	"complement, and work with NGS data'

setup(
	name='XYalign',
	version=__version__,
	description=package_description,
	author='Tim Webster',
	author_email='twebster17@gmail.com',
	url='https://github.com/WilsonSayresLab/XYalign',
	license='GNU General License version 3',
	long_description=open("README.md").read(),
	packages=["xyalign", "xyalign.test"],
	setup_requires=['pytest-runner'],
	tests_require=['pytest']
)
