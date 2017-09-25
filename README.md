# Renumbering
Python module for renumber molecules having same main structural framework with same number 

Getting started
------------
You need a setup.py file in the main directory.
You need pip which can be downloaded by  

```
sudo apt-get install python-pip 
```

pip can be later used to download the essential packages for python. Basically pip is a 'apt-get' or google store kind of thing for python libraries.

U need a __init__.py file in the main directory having all the codes. One is included and you can have a look.

To install the module, 
```
pip install -e .
```
on the main directory.

##Testing
Now install pytest for python 2.7 
In the test folder make a test file. One is already made with a assert true so that it passes now.
We need more tests.
```
sudo apt install python-logilab-common
sudo apt-get  install python-pytest
```

Finally to test the test:
```
python -m pytest test_orient.py
```

