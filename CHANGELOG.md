# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed 

### Removed 

### Fixed

## [0.2.0] - 2019-03-31
### Added
- Added `download_all_data()` function to `Universe` module to allow user to update
package data at runtime.
- Added solar flux and geomagnetic index data classes to `Universe` module.
- Added `EarthEnvironment` sub-package to store Earth-specific models and functionality
- Added `SpaceWeather` to `EarthEnvironment` module to store data releated to space weather and geomagnetic activity.
- Added `NRLMSISE00` module to `EarthEnvironment` package to implement higher-fidelity atmospheric density simulation.
- Added `Integrators` module to implement internal numerical integraiton in a 
manner which makes future improvements more easily accomplished.
- Added `Propagators` module to provide better dynamics integration tools.

### Changed
- Moved orbit dyanmics function `deriv_orbit_earth` out of `OrbitDynamics` module and into `Propagators` module. Also renamed function to `fderiv_earth_orbit`

### Removed
- Remove base `Simulation` module and associated funcitons including `propagate_orbit`

### Fixed
- Fixed scripts and functions used to download package data
- Fixed how many functions were specifying types and default arguments. There were
a number of instances of ineffectual type specification.
- Lots of documentation improvements and fixes.

## [0.1.1] - 2019-01-02
### Added

### Changed

### Removed

### Fixed
- Fixed `REQUIRE` to match `Project.toml` to address installation warning and 
package distribution issues.


## [0.1.0] - 2019-01-02
### Added
- Initial Release 

### Changed

### Removed

### Fixed