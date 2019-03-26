# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Added NRLMSISE00 Atmospheric model.
- Added solar flux and geomagnetic index data classes to `Universe` module.
- Added `download_all_data()` function to `Universe` module to allow user to update
package data at runtime.

### Changed

### Removed

### Fixed
- Fixed scripts and functions used to download package data
- Fixed how many functions were specifying types and default arguments. There were
a number of instances of ineffectual type specification.

## [0.1.0] - 2019-01-02
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