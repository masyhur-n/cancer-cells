# References

The citation keys match the `source` fields in `mgsim/parameters.py` and the
comments in `mgsim/model.py`. A value tagged **ASSUMPTION** in the code has no
direct measurement behind it. It is a modelling choice and is sampled over a
wide range in the Monte Carlo.

> Bibliographic details come from publisher and PubMed/PMC records. Check
> volume, page and DOI against the publisher before you cite any of them
> formally.

## Mg corrosion chemistry and in vivo degradation

| Key | Reference | Used for |
|---|---|---|
| Song1999 | Song G, Atrens A. Corrosion mechanisms of magnesium alloys. *Adv Eng Mater* 1999;1(1):11–33. | Reaction stoichiometry (1 Mg²⁺ : 2 OH⁻ : 1 H₂); Mg(OH)₂ pH ceiling |
| Song2007 | Song G. Control of biodegradation of biocompatable magnesium alloys. *Corros Sci* 2007;49(4):1696–1701. doi:10.1016/j.corsci.2007.01.001 | Tolerable H₂ ≈ 0.01 mL cm⁻² day⁻¹ (lower bound of `h2_tol`) |
| Staiger2006 | Staiger MP, Pietak AM, Huadmai J, Dias G. Magnesium and its alloys as orthopedic biomaterials: a review. *Biomaterials* 2006;27(9):1728–1734. doi:10.1016/j.biomaterials.2005.10.003 | Mg density and general properties |
| Witte2005 | Witte F, Kaese V, Haferkamp H, et al. In vivo corrosion of four magnesium alloys and the associated bone response. *Biomaterials* 2005;26(17):3557–3563. doi:10.1016/j.biomaterials.2004.09.049 | In vivo rates of AZ31 and WE43; bone response |
| Witte2010 | Witte F. The history of biodegradable magnesium implants: a review. *Acta Biomater* 2010;6(5):1680–1692. doi:10.1016/j.actbio.2010.02.028 | Clinical tolerance of gas (upper end of `h2_tol`) |
| Kraus2012 | Kraus T, Fischerauer SF, Hänzi AC, et al. Magnesium alloys for temporary implants in osteosynthesis: in vivo studies of their degradation and interaction with bone. *Acta Biomater* 2012;8(3):1230–1238. doi:10.1016/j.actbio.2011.11.008 | In vivo degradation of WE-type alloys |
| Kirkland2012 | Kirkland NT, Birbilis N, Staiger MP. Assessing the corrosion of biodegradable magnesium implants: a critical review of current methodologies and their limitations. *Acta Biomater* 2012;8(3):925–936. doi:10.1016/j.actbio.2011.11.014 | Gap between in vitro and in vivo rates (why the ranges are wide) |
| Zheng2014 | Zheng YF, Gu XN, Witte F. Biodegradable metals. *Mater Sci Eng R* 2014;77:1–34. doi:10.1016/j.mser.2014.01.001 | Alloy corrosion-rate ranges (HP-Mg, AZ31, ZK60, Mg-Zn-Ca) |
| Lamaka2018 | Lamaka SV, Gonzalez J, Mei D, et al. Local pH and its evolution near Mg alloy surfaces exposed to simulated body fluids. *Adv Mater Interfaces* 2018;5:1800169. doi:10.1002/admi.201800169 | Interfacial pH 9.9–10.5 without Ca, 7.8–8.5 with Ca; a Ca-P layer holds pH below about 8.5 (sets `ph_cap` for HA coatings) |

## Coatings (PEO / MAO, HA, Ca-P)

| Key | Reference | Used for |
|---|---|---|
| Narayanan2014 | Narayanan TSNS, Park IS, Lee MH. Strategies to improve the corrosion resistance of microarc oxidation (MAO) coated magnesium alloys for degradable implants: prospects and challenges. *Prog Mater Sci* 2014;60:1–71. doi:10.1016/j.pmatsci.2013.08.002 | PEO barrier factor, breakdown, sealing |
| Matykina2016 | Matykina E, Garcia I, Arrabal R, et al. Role of PEO coatings in long-term biodegradation of a Mg alloy. *Appl Surf Sci* 2016;389:810–823. [record](https://research.manchester.ac.uk/en/publications/role-of-peo-coatings-in-long-term-biodegradation-of-a-mg-alloy/) | About 35× lower H₂ for the first 2 weeks, steady state by week 6 (PEO `m0`, `t_break`) |
| Gu2011 | Gu XN, Li N, Zhou WR, et al. Corrosion resistance and surface biocompatibility of a microarc oxidation coating on a Mg–Ca alloy. *Acta Biomater* 2011;7(4):1880–1889. doi:10.1016/j.actbio.2010.11.034 | PEO on Mg-Ca |
| FattahAlhosseini2023 | Fattah-alhosseini A, et al. Review of HAp-containing PEO coatings on Mg. *J Magnes Alloys* 2023;11(9):2999–3011. [DOAJ](https://doaj.org/article/96a2ab2be0cb4f1ca47b32da83461dd7) | PEO+HA barrier and bioactivity |
| Hornberger2012 | Hornberger H, Virtanen S, Boccaccini AR. Biomedical coatings on magnesium alloys – a review. *Acta Biomater* 2012;8(7):2442–2455. doi:10.1016/j.actbio.2012.04.012 | HA / Ca-P coatings, osteoconduction bonus |
| Song2010 | Song Y, Zhang S, Li J, Zhao C, Zhang X. Electrodeposition of Ca–P coatings on biodegradable Mg alloy: in vitro biomineralization behavior. *Acta Biomater* 2010;6(5):1736–1742. doi:10.1016/j.actbio.2009.12.020 | HA-only (electrodeposited) coating |

## Mg degradation and tumour cells

| Key | Reference | Used for |
|---|---|---|
| Zhang2012 | Zhang Y, Ren L, Li M, Lin X, Zhao H, Yang K. Preliminary study on cytotoxic effect of biodegradation of magnesium on cancer cells. *J Mater Sci Technol* 2012;28(9):769–772. [record](https://www.jmst.org/EN/Y2012/V28/I9/769) | U2-OS cells: high alkalinity was cytotoxic, Mg²⁺ alone was not (`pH50_kill`) |
| Li2014 | Li M, Ren L, Li L, et al. Cytotoxic effect on osteosarcoma MG-63 cells by degradation of magnesium. *J Mater Sci Technol* 2014;30(9):888–893. [article](https://www.jmst.org/article/2014/1005-0302-30-9-888.html) | Bare Mg strongly cytotoxic, MAO-coated Mg much weaker |
| Burmester2014 | Burmester A, Luthringer B, Willumeit R, Feyerabend F. Comparison of the reaction of bone-derived cells to enhanced MgCl₂-salt concentrations. *Biomatter* 2014;4:e967616. [PMC4578555](https://pmc.ncbi.nlm.nih.gov/articles/PMC4578555) | MG63, SaOS2 and U2OS cell numbers fall above 10 mM Mg²⁺ (`IC50_mg`) |
| Globig2020 | Globig P, et al. Optimizing an osteosarcoma-fibroblast coculture model to study antitumoral activity of magnesium-based biomaterials. *Int J Mol Sci* 2020. [PMC7404313](https://pmc.ncbi.nlm.nih.gov/articles/PMC7404313/) | Relative tumour-cell reduction on Mg vs Ti |
| Globig2022 | Globig P, et al. Slow degrading Mg-based materials induce tumor cell dormancy on an osteosarcoma-fibroblast coculture model. *Bioact Mater* 2022. doi:10.1016/j.bioactmat.2021.12.031 | Cytostatic or dormancy response with no apoptosis at low rates; pH is the main driver (`pH50_static`, `s_max`) |
| Hakimi2016 | Hakimi O, Ventura Y, Goldman J, Vago R, Aghion E. Porous biodegradable EW62 medical implants resist tumor cell growth. *Mater Sci Eng C* 2016;61:516–525. doi:10.1016/j.msec.2015.12.043 | Faster-corroding and porous Mg is more cytotoxic (`area_factor`) |
| Qiao2020 | Qiao S, Wang Y, Zan R, et al. Biodegradable Mg implants suppress the growth of ovarian tumor. *ACS Biomater Sci Eng* 2020;6:1755–1763. | In vivo evidence of antitumour effect |
| Zan2021 | Zan R, et al. Biodegradable magnesium implants: a potential scaffold for bone tumor patients. *Sci China Mater* 2021. doi:10.1007/s40843-020-1509-2 | Concept of a dual-function Mg scaffold for bone tumours |
| Wu2016 | Unique antitumor property of the Mg-Ca-Sr alloys with addition of Zn. *Sci Rep* 2016. [PMC4764862](https://pmc.ncbi.nlm.nih.gov/articles/PMC4764862/) | Mg-Zn-Ca-type alloy antitumour effect |
| Review2024 | Mg alloys with antitumor and anticorrosion properties for orthopedic oncology: a review from mechanisms to application strategies. 2024. [PMC11026114](https://pmc.ncbi.nlm.nih.gov/articles/PMC11026114/) | Overview of mechanisms and of the corrosion-rate vs antitumour trade-off |

## Mg²⁺, pH and bone formation

| Key | Reference | Used for |
|---|---|---|
| Yoshizawa2014 | Yoshizawa S, Brown A, Barchowsky A, Sfeir C. Magnesium ion stimulation of bone marrow stromal cells enhances osteogenic activity, simulating the effect of magnesium alloy degradation. *Acta Biomater* 2014;10(6):2834–2842. doi:10.1016/j.actbio.2014.02.002 | Osteogenic stimulation by Mg²⁺ (`mg_stim_*`) |
| Zhang2016 | Zhang Y, Xu J, Ruan YC, et al. Implant-derived magnesium induces local neuronal production of CGRP to improve bone-fracture healing in rats. *Nat Med* 2016;22(10):1160–1169. doi:10.1038/nm.4162 | In vivo osteogenic benefit of Mg²⁺ |
| Wang2015 | Wang J, Witte F, Xi T, et al. Recommendation for modifying current cytotoxicity testing standards for biodegradable magnesium-based materials. *Acta Biomater* 2015;21:237–249. doi:10.1016/j.actbio.2015.04.011 | Toxicity of high-Mg extracts (`mg_tox_IC50`) |
| Galow2017 | Galow AM, Rebl A, Koczan D, Bonk SM, Baumann W, Gimsa J. Increased osteoblast viability at alkaline pH in vitro provides a new perspective on bone regeneration. *Biochem Biophys Rep* 2017;10:17–25. doi:10.1016/j.bbrep.2017.02.001 | MC3T3-E1 viability optimum near pH 8.4, mineralisation only at alkaline pH (`ph_opt`, `ph_amp`) |
| JahnenDechent2012 | Jahnen-Dechent W, Ketteler M. Magnesium basics. *Clin Kidney J* 2012;5(Suppl 1):i3–i14. doi:10.1093/ndtplus/sfr163 | Serum Mg 0.7–1.05 mM (baseline) |

## Tumour pH and Mg biology

| Key | Reference | Used for |
|---|---|---|
| Webb2011 | Webb BA, Chimenti M, Jacobson MP, Barber DL. Dysregulated pH: a perfect storm for cancer progression. *Nat Rev Cancer* 2011;11(9):671–677. doi:10.1038/nrc3110 | Tumour extracellular pH 6.7–7.1 (`dpH_acid`) |
| Swietach2014 | Swietach P, Vaughan-Jones RD, Harris AL, Hulikova A. The chemistry, physiology and pathology of pH in cancer. *Phil Trans R Soc B* 2014;369:20130099. doi:10.1098/rstb.2013.0099 | Bicarbonate buffering capacity (`beta`), tumour acidosis |
| GatenbyGawlinski1996 | Gatenby RA, Gawlinski ET. A reaction-diffusion model of cancer invasion. *Cancer Res* 1996;56(24):5745–5753. | Modelling framework: acid production by tumour cells, reaction–diffusion |
| WolfTrapani2008 | Wolf FI, Trapani V. Cell (patho)physiology of magnesium. *Clin Sci* 2008;114(1):27–35. doi:10.1042/CS20070129 | Mg²⁺ as a proliferation co-factor ("Mg-stimulated" phenotype) |
| CastiglioniMaier2011 | Castiglioni S, Maier JAM. Magnesium and cancer: a dangerous liason. *Magnes Res* 2011;24(3):S92–S100. doi:10.1684/mrh.2011.0285 | Mg²⁺ pro-tumour possibility |
