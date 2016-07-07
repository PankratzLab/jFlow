# Change Log
All notable changes to Genvisis will be documented in this file.

##[Unreleased]


##[0.1.4]
Reorganized mitoPipeline

##[0.1.3]
### Changed

- Default plotting types for mtDNA PC evaluation, now does individual plots for each variable. 

### Fixed

- PC correction iterator reporting one less PC than specified

## [0.1.2]
### Added
- Finalized Beta optimization in [MitoPipeLine](https://github.com/npankrat/Genvisis/commits/master/src/cnv/manage/MitoPipeline.java)
- Long format parsing is multi-threaded if more than one long format file is present

### Fixed
- Fixed issue where MitoPipeline was not downloading beta files


## [0.1.1]
### Added

- Possible fix to [issue](https://github.com/npankrat/Genvisis/issues/8), checking for NullPointerException when reading src files
- Also, addressed the [issue](https://github.com/npankrat/Genvisis/issues/8) where src files could fail import without a final warning.
- Automatic update check when GUI is launched

## 0.1.0 - 2016-05-20
### Added
- Initial Official versioning
### Changed 
- Versioning given by the build.xml file
### Removed
- an example of the removed heading

 
[Unreleased]: https://github.com/npankrat/Genvisis/compare/v0.1.3...HEAD
[0.1.4]: https://github.com/npankrat/Genvisis/compare/v0.1.2...v0.1.3
[0.1.3]: https://github.com/npankrat/Genvisis/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/npankrat/Genvisis/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/npankrat/Genvisis/compare/v0.1.0...v0.1.1

