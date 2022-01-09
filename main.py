from matplotlib import pyplot as plt
from pyopenms import *
import matplotlib.patches as mpatches
exp = MSExperiment()
#Load mzML file and store it in MSExperiment
MzMLFile().load("20190715_Redox_PISA_lysate_NADPH_set2-5.mzML",exp)
protein_ids = []
peptide_ids = []
SimpleSearchEngineAlgorithm().search("20190715_Redox_PISA_lysate_NADPH_set2-5.mzML","sequence.fasta", protein_ids, peptide_ids)
for peptide_id in peptide_ids:
  # Peptide identification values
  print(35*"=")
  print("Peptide ID m/z:", peptide_id.getMZ())
  print("Peptide ID rt:", peptide_id.getRT())
  print("Peptide scan index:", peptide_id.getMetaValue("scan_index"))
  print("Peptide ID score type:", peptide_id.getScoreType())

  # PeptideHits
  for hit in peptide_id.getHits():
    print(" - Peptide hit sequence:", hit.getSequence())
    mz = hit.getSequence().getMonoWeight(Residue.ResidueType.Full, hit.getCharge()) / hit.getCharge()
    print(" - Peptide hit monoisotopic m/z:", mz)
    print(" - Peptide ppm error:", abs(mz - peptide_id.getMZ())/mz *10**6 )
    print(" - Peptide hit score:", hit.getScore())
    score=mz/peptide_id.getMZ()
    print("comparsion score = ", score)
    a = str(hit.getSequence())
    tsg = TheoreticalSpectrumGenerator()
    spec1 = MSSpectrum()
    peptide = AASequence.fromString(a)
    p = Param()
    p.setValue("add_b_ion", "true")
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)
    tsg.getSpectrum(spec1, peptide, 1, 2)
    print("Spectrum 1 of", peptide, "has", spec1.size(), "peaks.")
    tsg1 = TheoreticalSpectrumGenerator()
    theo_spectrum = MSSpectrum()
    p = tsg1.getParameters()
    p.setValue("add_y_ions", "true")
    p.setValue("add_b_ions", "true")
    p.setValue("add_metainfo", "true")
    tsg1.setParameters(p)
    peptide = AASequence.fromString(a)
    tsg1.getSpectrum(theo_spectrum, peptide, 1, 2)
    experiment_spectrum = exp.getSpectrum(peptide_id.getMetaValue("scan_index"))
    alignment = []
    spa = SpectrumAlignment()
    p = spa.getParameters()
    # use 0.5 Da tolerance (Note: for high-resolution data we could also use ppm by setting the is_relative_tolerance value to true)
    p.setValue("tolerance", 0.5)
    p.setValue("is_relative_tolerance", "false")
    spa.setParameters(p)
    # align both spectra
    spa.getSpectrumAlignment(alignment, theo_spectrum, experiment_spectrum)
    print("Number of matched peaks: " + str(len(alignment)))
    print("ion\ttheo. m/z\tobserved m/z")

    for theo_idx, obs_idx in alignment:
        ion_name = theo_spectrum.getStringDataArrays()[0][theo_idx].decode()
        ion_charge = theo_spectrum.getIntegerDataArrays()[0][theo_idx]
        print(ion_name + "\t" + str(ion_charge) + "\t"
              + str(theo_spectrum[theo_idx].getMZ())
              + "\t" + str(experiment_spectrum[obs_idx].getMZ()))

    theo_mz, theo_int, obs_mz, obs_int = [], [], [], []
    for theo_idx, obs_idx in alignment:
        theo_mz.append(theo_spectrum[theo_idx].getMZ())
        theo_int.append(theo_spectrum[theo_idx].getIntensity())
        obs_mz.append(experiment_spectrum[obs_idx].getMZ())
        obs_int.append(experiment_spectrum[obs_idx].getIntensity())
    plt.subplot(2, 1, 1)
    # make range form 300 to 1000 to see the simialrity

    for mz, intensity in zip(*experiment_spectrum.get_peaks()):
        if mz >= min(obs_mz) and mz <= max(obs_mz):
            plt.plot([obs_mz, obs_mz], [0, intensity], color='red')
            plt.plot([mz, mz], [0, intensity], color='black')
    b=str(peptide_id.getMetaValue("scan_index"))
    plt.title(a+" "+b)
    plt.subplot(2, 1, 2)
    color_bar = ['red', 'blue']
    count = 0
    for mz, i in zip(*theo_spectrum.get_peaks()):
        if mz >= min(theo_mz) and mz <= max(theo_mz):
            if 'y' in theo_spectrum.getStringDataArrays()[0][count].decode():
                plt.plot([mz, mz], [0, i], snap=False, color=color_bar[0])
                count = count + 1
            else:
                plt.plot([mz, mz], [0, i], snap=False, color=color_bar[1])
                count = count + 1

    red_patch = mpatches.Patch(color='red', label='y ions')
    blue_patch = mpatches.Patch(color='blue', label='b ions')
    plt.legend(handles=[red_patch, blue_patch], loc='lower right')
    plt.ylim(bottom=0)
    plt.show()
