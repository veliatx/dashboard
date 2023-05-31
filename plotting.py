import matplotlib.pyplot as plt

def plot_structure_plddt(plddt):
    fig, ax = plt.subplots(figsize=(4, 3))
    fs=8
    n_aa = len(plddt)
    ax.plot(range(n_aa), plddt)
    ax.set_title("ESMFold Prediction Confidence", fontsize=fs)
    ax.set_xlabel("Amino Acid", fontsize=fs)
    ax.set_ylabel("plDDT", fontsize=fs)
    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), fontsize=fs)
    ax.set_yticks(ax.get_yticks(), ax.get_yticklabels(), fontsize=fs)
    ax.set_xlim(0, n_aa)
    ax.set_ylim([0, 100])
    ax.hlines(50, 0, n_aa, linestyle='--', color='r', alpha=0.5, label='Low Confidence')
    ax.hlines(90, 0, n_aa, linestyle='--', color='g', alpha=0.5, label='High Confidence')
    return fig, ax
