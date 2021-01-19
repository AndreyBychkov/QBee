import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="qbee-AndreyBychkov",
    version="0.5.0",
    author="Andrey Bychkov",
    author_email="abychkov@edu.hse.ru",
    description="Package for transforming ODEs system to quadratic form.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/AndreyBychkov/QBee",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)