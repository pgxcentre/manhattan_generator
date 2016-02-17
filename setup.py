#!/usr/bin/env python

# How to build source distribution
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip


from setuptools import setup


VERSION="1.7.1"


def setup_package():
    setup(
        name="manhattan_generator",
        version=VERSION,
        description="Creation of beautiful Manhattan plots",
        author="Louis-Philippe Lemieux Perreault",
        author_email="louis-philippe.lemieux.perreault@statgen.org",
        url="https://github.com/pgxcentre/manhattan_generator",
        license="CC BY-NC 4.0",
        entry_points={
            "console_scripts": [
                "manhattan_generator=manhattan_generator:safe_main",
            ],
        },
        py_modules=["manhattan_generator"],
        install_requires=["matplotlib >=1.3.1", "numpy >= 1.8.0",
                          "pandas >= 0.17.0"],
        classifiers=[
            "Operating System :: Linux",
            "Programming Language :: Python",
            "Programming Language :: Python :: 2.7",
        ],
    )

    return


if __name__ == "__main__":
    setup_package()
