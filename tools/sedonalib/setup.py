from setuptools import setup, find_packages

setup(name='sedonalib',
      version='1.0.0',
      description='Sedona python library',
      author='sedonacrew',
      author_email='email',
      license='whatever, e.g. BSD',
      packages=find_packages(),
      package_data={"sedonalib": ["models/*", "params/lightcurve/*",
                                  "plt/*", "spectrum/*"]},
      install_requires=['numpy', 'scipy', 'h5py','matplotlib==2.1.0'],
      zip_safe=False)
