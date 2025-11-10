# TEFs for selected groups of compounds

A list of named vectors which provide Toxic Equivalency Factors (TEFs)
for calculating Toxic EQuivalent (TEQ) sums.

## Usage

``` r
info_TEF$DFP_2005
info_TEF$DFP_2022
```

## Format

A list of named vectors:

- `DFP_2005` gives the 2005 World Health Organisation TEFs for dioxins,
  furans and dioxin-like (planar) polychlorinated biphenyls ([Van den
  Berg et al. 2006](https://doi.org/10.1093/toxsci/kfl055))

- `DFP_2022` gives the 2022 World Health Organisation TEFs for dioxins,
  furans and dioxin-like (planar) polychlorinated biphenyls ([DeVito et
  al. 2024](https://doi.org/10.1016/j.yrtph.2023.105525))

- `DFP_HOLAS3` is a superseded version of `DFP_2005` that was used in
  the HELCOM HOLAS3 assessment. It excludes CB114, CB123 and CB189 and
  uses code CDFO instead of OCDF. It is included for reproducibility,
  but might be removed in the future.

- `DFP_CEMP` is a superseded version of `DFP_2005` that was used in the
  OSPAR CEMP assessments up to and including 2024. It excludes CB114,
  CB123 and CB189, uses code CDFO instead of OCDF, and has a TEF for
  CDFO of 0.00003 rather than 0.0003. It is included for
  reproducibility, but might be removed in the future.

## References

DeVito M, Bokkers B, van Duursen MBM, van Ede K, Feeley M, Gáspár EAF,
Haws L, Kennedy S, Peterson RE, Hoogenboom R, Nohara K, Petersen K,
Rider C, Rose M, Safe S, Schrenk D, Wheeler MW, Wikoff DS, Zhao B, van
den Berg M, 2024. The 2022 world health organization reevaluation of
human and mammalian toxic equivalency factors for polychlorinated
dioxins, dibenzofurans and biphenyls. Regulatory Toxicology and
Pharmacology 146; <https://doi.org/10.1016/j.yrtph.2023.105525>

Van den Berg M, Birnbaum LS, Denison M, De Vito M, Farland W, Feeley M,
Fiedler H, Hakansson H, Hanberg A, Haws L, Rose M, Safe S, Schrenk D,
Tohyama C, Tritscher A, Tuomisto J, Tysklind M, Walker N, Peterson RE,
2006. The 2005 World Health Organization reevaluation of human and
mammalian toxic equivalency factors for dioxins and dioxin-like
compounds. Toxicological Sciences 93 223–241;
<https://doi.org/10.1093/toxsci/kfl055>
