var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#SatelliteDynamics.jl-Documentation-1",
    "page": "Home",
    "title": "SatelliteDynamics.jl Documentation",
    "category": "section",
    "text": "Welcome to the SatelliteDynamics.jl documentation! The package is meant to address the needs of the satellite operator, academic researcher, and public enthusiast  communities. In particular, it is designed and built to fill the following needs:Open-source, MIT-lencesed software to remove barriers to entry of doing high-fidelity satellite modeling.\nHigh-quality, validated, tested library for orbit and attitude dynamics modeling.\nEasily acceible API design to make implementation of simulation and analysis intuitive."
},

{
    "location": "#Supporting-and-Citing-1",
    "page": "Home",
    "title": "Supporting and Citing",
    "category": "section",
    "text": "The software in this ecosystem was developed as part of academic research. If you would like to help support it, please star the repository as such metrics on usage help guage interest and secure funding. If you use the software as part of your research, teaching, or other activities, we would be grateful if you could cite our work."
},

{
    "location": "#Getting-Started:-Installation-And-First-Steps-1",
    "page": "Home",
    "title": "Getting Started: Installation And First Steps",
    "category": "section",
    "text": "To install the package, use the following command inside the Julia REPL:Pkg.add(\"SatelliteDynamics\")To load the package, use the command:using SatelliteDynamics"
},

{
    "location": "#Package-Structure-1",
    "page": "Home",
    "title": "Package Structure",
    "category": "section",
    "text": "The package is divided into a number of submodules each designed to provide a single, well-defined set of functions. The details on Pages = [\n    \"modules/constants.md\",\n    \"modules/universe.md\",\n]\nDepth = 2"
},

{
    "location": "#Examples-1",
    "page": "Home",
    "title": "Examples",
    "category": "section",
    "text": "The best way to learn how to use any software is to ``@contents Pages = [     \"tutorials/universeexample.md\",     \"tutorials/timeexample.md\",     \"tutorials/epoch_example.md\", ] Depth = 2\n## Index\n@index ```"
},

{
    "location": "modules/constants/#",
    "page": "Constants",
    "title": "Constants",
    "category": "page",
    "text": ""
},

{
    "location": "modules/constants/#Constants-1",
    "page": "Constants",
    "title": "Constants",
    "category": "section",
    "text": "The constants module constains common mathematical, physical, astronomical constants."
},

{
    "location": "modules/universe/#",
    "page": "Universe",
    "title": "Universe",
    "category": "page",
    "text": ""
},

{
    "location": "modules/universe/#Universe-1",
    "page": "Universe",
    "title": "Universe",
    "category": "section",
    "text": "The Universe module defines simulation-specific data files which are constants  of most simulations. In particular it provides data structures for storing and accessing Earth orientation parameters and spherical harmonic gravity field  models.The module defines the global variables EOP and GRAVITY_MODEL which are loaded at runtime and .EOP defaults to use the rapid Earth orientation data file finals.all (IAU 2000) distributed by the IERS. The module also supports IERS C04 product files.GRAVITY_MODEL defaults to use the EGM2008 spherical harmonic gravity model,  truncated to order and degree 90."
},

{
    "location": "tutorials/universe_example/#",
    "page": "Earth Orientation Parameters",
    "title": "Earth Orientation Parameters",
    "category": "page",
    "text": ""
},

{
    "location": "tutorials/universe_example/#Earth-Orientation-Parameters-1",
    "page": "Earth Orientation Parameters",
    "title": "Earth Orientation Parameters",
    "category": "section",
    "text": ""
},

{
    "location": "tutorials/epoch_example/#",
    "page": "Time Epoch",
    "title": "Time Epoch",
    "category": "page",
    "text": ""
},

{
    "location": "tutorials/epoch_example/#Time-Epoch-1",
    "page": "Time Epoch",
    "title": "Time Epoch",
    "category": "section",
    "text": ""
},

]}
