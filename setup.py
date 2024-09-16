from setuptools import setup, Extension

# Define the extension module
symnmf_module = Extension('symnmfmodule',
                          sources=['symnmfmodule.c', 'symnmf.c'],
                          include_dirs=[],
                          libraries=[],  
                          library_dirs=[],  
                          )

# Setup function
setup(
    name='SymNMF',
    version='0.1',
    description='Python Package with C extension for Symmetric NMF',
    ext_modules=[symnmf_module]
)
