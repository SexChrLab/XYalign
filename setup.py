from setuptools import setup

# Get version from xyalign/_version.py

version_py = "xyalign/_version.py"
exec(open(version_py).read())

package_description = 'Command line tools and python library to infer ploidy,'\
	' correct for sex chromosome complement, and work with NGS data'

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
	tests_require=['pytest'],
	install_requires=[
		'matplotlib', 'numpy', 'pandas', 'pybedtools', 'pysam', 'scipy'],
	py_modules=[
		'xyalign.xyalign', 'xyalign.scripts.explore_thresholds',
		'xyalign.scripts.plot_count_stats',
		'xyalign.scripts.plot_window_differences',
		'xyalign.scripts.compare_vcfs'],
	entry_points={
		"console_scripts": [
			'xyalign = xyalign.xyalign:main',
			'explore_thresholds = xyalign.scripts.explore_thresholds:main',
			'plot_count_stats = xyalign.scripts.plot_count_stats:main',
			'plot_window_differences = xyalign.scripts.plot_window_differences:main',
			'compare_vcfs = xyalign.scripts.compare_vcfs:main']}
)
