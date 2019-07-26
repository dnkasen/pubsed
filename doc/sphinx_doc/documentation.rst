================================
Adding Documentation
================================

This documentation was constructed from .rst files in the
**docs/sphinx_doc/** directory

For basics of using the reStructuredTex (.rst) markup language, see

http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

To compile the documentation into html or latex you must install Sphinx from

http://www.sphinx-doc.org/en/master/usage/installation.html

e.g., on a mac you can try to install using::

  brew install sphinx-doc

Once installed, make the webpage version from the sphinx_doc/ directory using::

  make html

The resulting webpage files appear in the _build/html/ directory.
You can view them using::

  open _build/html/index.html

To make a latex file use::

  make latex

which makes a complete latex file in _build/latex/. To make the latex file and
compile it into a pdf using pdflatex you can just use::

  make latexpdf
