#!/usr/bin/env bash
rm -rf build dist treeswift.*
python3 setup.py sdist
python3 setup.py bdist_wheel --universal
twine upload dist/*
