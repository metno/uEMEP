---
project: uEMEP
summary: Air quality dispersion model for high resolution downscaling of EMEP-MSC-W
author: Norwegian Meteorological Institute
website: https://www.met.no/en
twitter: https://x.com/meteorologene
github: https://github.com/metno
project_github: https://github.com/metno/uEMEP
src_dir: ./src
output_dir: ./doc
source: true
incl_src: true
graph: true
search: true
display: public
         private
         protected
preprocess: false
page_dir: ./docs
sort: alpha
version: 7.0.0
parallel: 8
---

The uEMEP (urban EMEP) model was created to extend the EMEP MSC-W model's application to near-street-level air quality modeling. It uses a downscaling approach based on classical Gaussian plume modeling, which is combined with the physical parameterization and emission data from the EMEP MSC-W model.
