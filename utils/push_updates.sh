#!/usr/bin/env bash

if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    # Install required packages for coverage and documentation
    julia --project -e 'import Pkg; Pkg.add("Coverage");'
    julia --project -e 'import Pkg; Pkg.add("Documenter");'
    julia --project -e 'import Pkg; Pkg.add("Plots");'
    julia --project -e 'import PKG; Pkg.add("LaTeXStrings")'

    # Submit test coverage report
    julia --project -e 'using Coverage; Coveralls.submit(Coveralls.process_folder())'

    # Build and deploy documentation
    julia --project ./docs/make.jl
fi