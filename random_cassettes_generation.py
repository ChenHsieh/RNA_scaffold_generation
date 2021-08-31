# -*- coding: utf-8 -*-
"""random_cassettes_generation.py

# sequence element input
"""

import seaborn as sns
from collections import Counter
import random
import pandas as pd
import RNA
run_name = "Chen-testing-1"
# update gRNA every run, the total generated length would be the same as the gRNA number provided
gRNAs = {
    "g1": "AACGATGTGTTGGTTACACAAGG",
    "g2": "CTTCATACTGTTGGGCCAAAGGG",
    "g3": "CTTTCTACTCAAATCCAGTTTGG",
    "g4": "AATGGTGCCTTGTAGAAGGTGGG",
    "g5": "TTACTATAACATATAGATAGAGG",
    "g6": "ACGTGGTTTCAGATTACACGTGG"
}

# default sequence element collection
promoters = {"MtU6": "ATGCCTATCTTATATGATCAATGAGGCATTTAATTGGGTGCATATGATGGTGAAAAAAGGTGCAGCTCCTGGCTTGGGAATGATGACTCATGTGGAATTTGGTCTTAAATTTATCACATCCTTTTGGGATGTGATGATTGTATCACTTGTTCATTTTGCAAAGACAAGGTGCACTGCTACAAACTTTGGTTTAATCTGAAATAAAACAAAACTCACTGAGAGGAAGATGCATCCCAGTAGGTGAAAGTCGAGAAGGATTTGCATGTTACTATTACACTTGCTTTTTAGTCCCACATCGTCTGAAACATAAAATATTTCAGCGTTTAAATACTTCAAGCGAACCAGTAGGCTTG", "AtU6.26P205": "CATCTTCATTCTTAAGATATGAAGATAATCTTCAAAAGGCCCCTGGGAATCTGAAAGAAGAGAAGCAGGCCCATTTATATGGGAAAGAACAATAGTATTTCTTATATAGGCCCATTTAAGTTGAAAACAATCTTCAAAAGTCCCACATCGCTTAGATAAGAAAACGAAGCTGAGTTTATATACAGCTAGAGTCGAAGTAGTGATTG", "AtU6.26(181)": "ATAATCTTCAAAAGGCCCCTGGGAATCTGAAAGAAGAGAAGCAGGCCCATTTATATGGGAAAGAACAATAGTATTTCTTATATAGGCCCATTTAAGTTGAAAACAATCTTCAAAAGTCCCACATCGCTTAGATAAGAAAACGAAGCTGAGTTTATATACAGCTAGAGTCGAAGTAGTGATTG", "AtU6-26(79)": "GGAGTGATCAAAAGTCCCACATCGATCAGGTGATATATAGCAGCTTAGTTTATATAATGATAGAGTCGACATAGCGATTG", "AtU6.1": "AGAAATCTCAAAATTCCGGCAGAACAATTTTGAATCTCGATCCGTAGAAACGAGACGGTCATTGTTTTAGTTCCACCACGATTATATTTGAAATTTACGTGAGTGTGAGTGAGACTTGCATAAGAAAATAAAATCTTTAGTTGGGAAAAAATTCAATAATATAAATGGGCTTGAGAAGGAAGCGAGGGATAGGCCTTTTTCTAAAATAGGCCCATTTAAGCTATTAACAATCTTCAAAAGTACCACAGCGCTTAGGTAAAGAAAGCAGCTGAGTTTATATATGGTTAGAGACGAAGTAGTGATTG", "AtU3d(old)": "ATAAGCTTATGATTTCTTTTTTCTTACGAATTTTGCGTCCCACATCGGTAAGCGAGTGAAGAAATAACTGCTTTATATATGGCTACAAAGCACCATTGGTCA", "AtU3b": "TTATGTAATTTAATGGGCTATCGTCCATATATTCACTAATACCCATGCCCAGTACCCATGTATGCGTTTCATATAAGCTCCTAATTTCTCCCACATCGCTCAAATCTAAACAAATCTTGTTGTATATATAACACTGAGGGAGCACCATTGGTCA", "VvU3.1": "AGTACTTTCATAGGAATAGGTTTCAAATGAACCTTTGTGATACACTTCGATCCTGACTCTCTCTAAAAGCAAAAATCATTAATATTTTTTCAATTATATTTTAATTTTTACAAATACTAACAAATATTAATAATTCTAATATCTTTCTTGTTTTGAAAAATAAAAGAAAAATAATAATGTTTGGTATATTCCGTATATATTATTTTAAATGAATCCAGAAGTTTCCAAGAATTTCACTGGCAATCAATCGTGCATCAGCTGTCAATCGTTGTTCCCAGGAAGGCTCATTGGAAGTCTATAACCAATGAGAACACGCGGTGACTAGCCGTCCCACATCGAAAATGCAGGAAACATTTAATAACTATATAACAAAGGATAGGAGATTCACATGCCA",
             "VvU3.2": "TATTTTTAGTTTTGTTAAAAATAATTATTTTTTAGAAAAGTATTTATAATAAAAATATTATTAAAAAATATTTTACGAATTAAAAGTACTTTCATAGGAATAGAAAAAAATCATTTACCTTAATTTTTAAAAAACTTTTACTATGAGATAAACATGATTTTCTCAAAAAAATTTAATCGCAGAATCTTAAAAATAGGATTCAAGTGACAGATCCGTTACTTTACACGTACAACTATTATCTTCTTTAACGATATAAATAAAAAAATGTTTATCAAAATAATTTTGCAAATTAAAATATCATACATCGGCCGTTAAAATTATTTTAAATGAATCCAGAAGTTTGCAAGAATTTCATTTGGCAATCTGTCGTGCATCAGCTGTCAATCGTTGTTACCAAGAAGTTCATTGCAAATATATGGCCAATGAGAACACGCGGTGACTAGCCGTCCCACATCGGAAATGCAGGAAACATTTAATGACTATATAACAAAGGATAGGAGATTCACATGCCA", "VvU6.1": "TTGCCTCTGGAAAATCCCCCTATAATTTTTTCAGTTTCTTTCAATTACCAGCAAAATTAGGGCTGAGAATGGAATTCCCAAACCCTAAAACAAAAAGAATTCCCAAATCTGAAAACAGTCTAGGATTCGGTTTTATAAATTACAAAAAGAACAACCATAATAAAGATGAAAAGCAACGAAGGAAGAAAAAGGAAATGAGAGAGGAGAGAGGAGACGCAGAGAGAGCAGCAGTCTCACCTTCTGGTGGGGAACACCAAGGACGAACATGCCAATTCTAAATTCAACCCAAATGAGTTGTGGTGACGGGCCGTGGGCTCGATGCCCAGACCAAGCGAAACGACGTCGTTCCCAAACGACTCTTCCCACATCGACTGCTCATAGACGAAATTGAGCTTTTATATATCAGGAGCAAACGCTTAGAGCTTG", "VvU6.2": "TTCCCTAATCATCATGTCTCTGCAATCTCTATCTAAGAGCTCCGCAAGAGATGATATTCAAAATCTTTTTGCCATGGTTTTGCTTCAAGTTGGTTTTTGATCAGCAGTCAATGGATTTTAAGCTACCACCCTCGGTATCCTACATAAGAAATCCAATACAAAAGTGGATTTTTGCAGTGCTGGTAGTTTCTTGAATTTAAGTTATTTAATTAGACTTATGATAAACACTAGCCCCTGGAAAATTCACCTACAATTTTTTTCAGTTTCTTTCAATTACCAGCAAAATTAGGGCTGAGAATGCAATTCCCAAACCCTAAGAACAAAAAGAATGCCCAAATCTGAAAGCATTATATGACAACAACCATAATAAAGATGAAAAGCAACGAAGAAAAAGTAATCTCACCTTCTGGTGGGGAACACCAAGGACGAACATGCCAATTCTAAATTCAACCCAAATGAGTTGTGGTGACGGGCCGTGGGCTCGATGCCCAGACCAAGCGAAACGACGTCGTTCCCAAACGACTCTTCCCACATCGACTGCGTATAGACTAAATTCACCTTTTTTATATCAGGAGCAAACGCTTAGAGCTTG"}
# promoters = {
#     "MtU6": "ATGCCTATCTTATATGATCAATGAGGCATTTAATTGGGTGCATATGATGGTGAAAAAAGGTGCAGCTCCTGGCTTGGGAATGATGACTCATGTGGAATTTGGTCTTAAATTTATCACATCCTTTTGGGATGTGATGATTGTATCACTTGTTCATTTTGCAAAGACAAGGTGCACTGCTACAAACTTTGGTTTAATCTGAAATAAAACAAAACTCACTGAGAGGAAGATGCATCCCAGTAGGTGAAAGTCGAGAAGGATTTGCATGTTACTATTACACTTGCTTTTTAGTCCCACATCGTCTGAAACATAAAATATTTCAGCGTTTAAATACTTCAAGCGAACCAGTAGGCTT",
#     "AtU3d": "ATAAGCTTATGATTTCTTTTTTCTTACGAATTTTGCGTCCCACATCGGTAAGCGAGTGAAGAAATAACTGCTTTATATATGGCTACAAAGCACCATTGGTC",
#     "AtU6.26P205": "CATCTTCATTCTTAAGATATGAAGATAATCTTCAAAAGGCCCCTGGGAATCTGAAAGAAGAGAAGCAGGCCCATTTATATGGGAAAGAACAATAGTATTTCTTATATAGGCCCATTTAAGTTGAAAACAATCTTCAAAAGTCCCACATCGCTTAGATAAGAAAACGAAGCTGAGTTTATATACAGCTAGAGTCGAAGTAGTGATT",
#     "AtU6-26(79)": "GGAGTGATCAAAAGTCCCACATCGATCAGGTGATATATAGCAGCTTAGTTTATATAATGATAGAGTCGACATAGCGATT",
#     "AtU6.1": "AGAAATCTCAAAATTCCGGCAGAACAATTTTGAATCTCGATCCGTAGAAACGAGACGGTCATTGTTTTAGTTCCACCACGATTATATTTGAAATTTACGTGAGTGTGAGTGAGACTTGCATAAGAAAATAAAATCTTTAGTTGGGAAAAAATTCAATAATATAAATGGGCTTGAGAAGGAAGCGAGGGATAGGCCTTTTTCTAAAATAGGCCCATTTAAGCTATTAACAATCTTCAAAAGTACCACAGCGCTTAGGTAAAGAAAGCAGCTGAGTTTATATATGGTTAGAGACGAAGTAGTGATT"
# }
linkers = {
    "linker1": "GCTACCTCAGCATAGTCTACAGC",
    "linker2": "CCACATGGGGACACGTGGCGCT",
    "linker3": "GAGCTGTCTGCATGTGTGTCAGC",
    "linker4": "GGACTGCGATTATGGAGCGTGC"
}
scaffolds = {
    "SF (standard)": "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT",
    "SF1": "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT",
    "SF2": "GTTTAAGAGCTATGCTGAAAAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTAaCAACggGAAAccGTGGCACCGAGTCGGTGCTTTTTTT",
    "SF3": "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT",
    "SF4": "GTTTGAGAGCTATGCTGGAAACAGCATAGCAAGTTCAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTCGGTTTTTTT",
    "SF5 (SF4558)": "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCACGCCGAAAGGCGGGCACCGAGTCGGTGCCTTTTTTT",
    "SF5b (error)": "GTTTTAGAGCTAGAAACGGCAAGTTAAAATAAGGCTAGTCCGTTATCACGCCGAAAGGCGGGCGCCGAGTCGGCGCC",
    "SF5a (4558a)": "GTTTCAGAGCTAGAAATAGCAAGTTGAAATAAGGCTAGTCCGTTATCACGCCGAAAGGCGGGCACCGAGTCGGTGCCTTTTTTT",
    "SF4561": "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCACCTTGGGCCGAAAGGCCCAAGGGGCGCCGAGTCGGCGCCTTTTTTT",
    "SF1a": "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTaATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT",
    "SF4a": "GTTTGAGAGCTATGCTGGAAACAGCATAGCAAGTTCAAATAAGGCTAGTCCGTaATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTCGGTTTTTTT"
}

"""# parameters"""

# CRISPR system: No need for the 5’-G if switching to Pol-II single transcriptional unit strategy
# DNA element arrangement

# some examples
# MtU6-g5-SF2-linker2-AtU3d-g6-SF4558-linker3-AtU6.26(79)-g4-SF4a-linker1
# linker1--AtU3d-g1-SF4558a -linker3-AtU6.1-g2-SF2-linker4AtU6.26P205-g3-SF

# hypothetical product
# MtU6-g5-SF2-linker2-AtU3d-g6-SF4558a-linker3-AtU6.26(79)-g4-SF4a-linker1-AtU6.1-g3-SF2-linker4-AtU3d-g1-SF1 -linker3-AtU6.26P205-g2-SF

element_arrangement = {
    "promoter": promoters,
    "gRNA": gRNAs,
    "scaffold": scaffolds,
    "linker": linkers
}

"""# functions"""

# function for calculating kmer


def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1
    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)
    return kmers

# get structure and energy


def get_structure(seq):
    fc = RNA.fold_compound(seq)
    # compute MFE and MFE structure
    (mfe_struct, mfe) = fc.mfe()
    return (mfe_struct, mfe)
# some script copied from forgi to deal with file reading

# swap nucleotide at stem structure


def swapBp(seq, index):
    comp = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }
    # print(index)
    print(seq[0:index]+'\''+seq[index]+'\''+seq[index+1:])
    print(seq[0:index]+'\''+comp[seq[index]]+'\''+seq[index+1:])
    neo_seq = seq[0:index]+comp[seq[index]]+seq[index+1:]
    return neo_seq


class OutFiletype:
    def __init__(self, write_fun, extention_fun, rna_type):
        self.convert = write_fun
        self.get_extention = extention_fun
        self.rna_type = rna_type


def to_bg_or_cg_string(cg):
    if not hasattr(cg, "coords") or not cg.coords.is_filled:
        return cg.to_bg_string()
    else:
        return cg.to_cg_string()


def cg_or_bg_extention(cg):
    if not hasattr(cg, "coords") or not cg.coords.is_filled:
        return ".bg"
    else:
        return ".cg"


FILETYPES = {
    "forgi": OutFiletype(to_bg_or_cg_string, cg_or_bg_extention, "any")
}


# save the prediction result to file
def save_structure_fx(scaffold_id, seq):
    (mfe_struct, mfe) = get_structure(seq)
    fx = f""">{scaffold_id}_{mfe}
    {seq}
    {mfe_struct}"""
    file_name = f"{scaffold_id}.fx"
    with open(file_name, 'w') as fx_file:
        fx_file.write(fx)
    return file_name


def get_structure_definition(scaffold_id, seq):
    file_name = save_structure_fx(scaffold_id, seq)
    cg = forgi.load_rna(file_name, allow_many=False)
    target_str = FILETYPES["forgi"].convert(cg)
    # fig = plt.figure()
    # fvm.plot_rna(cg, text_kwargs={"fontweight": "black"}, lighten=0.1,
    #  backbone_kwargs={"linewidth": 3})
    # plt.show()
    # plt.savefig(f"{file_name}.png")
    # plt.close(fig)
    # print(target_str)
    define_structure = [line.split(' ') for line in target_str.split(
        '\n') if line.split(' ')[0] == 'define']
    return define_structure


def random_swap_stem(define_structure, seq, times):
    for define in define_structure:
        if define[1].startswith('s'):
            print(define)
            forward = seq[int(define[2])-1:int(define[3])]
            reverse = seq[int(define[4])-1:int(define[5])]
            print(forward)
            print(reverse)
            reverse = reverse[::-1]
            for i in range(times):
                swapping_offset_index = random.randint(0, len(forward))
                # print(swapping_offset_index)
                forward_swapping_bp = int(define[2])-1+swapping_offset_index
                reverse_swapping_bp = int(define[5])-swapping_offset_index-1
                # print(seq)
                seq = swapBp(seq, forward_swapping_bp)
                seq = swapBp(seq, reverse_swapping_bp)
                # print(seq)
    return seq


def check_dup_kmer(kmers):
    dup_kmer_dict = {}
    for kmer, count in Counter(kmers).items():
        if count > 1:
            dup_kmer_dict[kmer] = count
    return dup_kmer_dict


def draw_random_element(element_dict):
    entry_list = list(element_dict.items())
    (sign, sequence) = random.choice(entry_list)
    return (sign, sequence)


def generate_sub_cassette(gRNA_sign, gRNA, mode):
    # initiate variable
    sub_cassette_sign = []
    sub_cassette = ""
    # go through each element of cassette
    for element_type, element_dict in element_arrangement.items():
        # deal with gRNA, the only not random element
        if element_type == "gRNA":
            sub_cassette_sign.append(gRNA_sign)
            sub_cassette = sub_cassette + gRNA
            continue
        # fill up the cassette elements
        (sign, sequence) = draw_random_element(element_dict)
        sub_cassette_sign.append(sign)
        sub_cassette = sub_cassette + sequence
    return sub_cassette_sign, sub_cassette


def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1
    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)
    return kmers


def check_dup_kmer(kmers):
    dup_kmer_dict = {}
    for kmer, count in Counter(kmers).items():
        if count > 1:
            dup_kmer_dict[kmer] = count
    return dup_kmer_dict


def condition_check(hypothetical_product):

    # LENGTH 0.3kb - 1.8kb
    inRange = len(hypothetical_product) in range(300, 1800)

    # print(len(hypothetical_product))
    # 19 mer check
    kmers = build_kmers(hypothetical_product, 19)
    dup_kmers = check_dup_kmer(kmers)
    noKmer = not dup_kmers
    isPass = all([inRange, noKmer])

    return isPass


"""# generation"""

# init
gRNA_sign_list = list(gRNAs.keys())
gRNA_num = len(gRNA_sign_list)
product_number = 1000
hypothetical_products = {}
simulation_count = 0
product_count = 0
simulating = True
halt = 1000000
fragment_number = 3

fragment_element_length = len(gRNAs)/fragment_number*4-1

serial_number_list = []
part_list = []
content_sign_list = []
seq_list = []
length_list = []
mfe_list = []
mfe_struct_list = []

product_dict = {}

product_start = "CAAGCGAACCAGTAGGCTTG"
product_end = "GTTTTAGAGCTAGAAATAGC"

# start simulation
while simulating:
    # for controling the fixed product start
    starting_fragment = True
    simulation_count += 1

    hypothetical_product_sign = []
    hypothetical_product = ""

    synthesis_fragments_sign = []
    synthesis_fragments = []

    current_synthesis_fragment_sign = []
    current_synthesis_fragment = ""

    # the structure of the synthesis fragments

    # `fixed U6 promoter 3' 20bp`-gRNA-scaffold-linker promoter-gRNA-scaffold-linker
    #           linker + promoter-gRNA-scaffold-linker promoter-gRNA-scaffold-linker
    #           linker + promoter-gRNA-scaffold-linker promoter-gRNA-scaffold-linker
    #           ...
    #           linker + promoter-gRNA-scaffold-linker promoter-gRNA-`fixed SF 5' 20bp`

    # go through each gRNA randomly
    gRNA_count = 0
    fragment_count = 0
    random.shuffle(gRNA_sign_list)
    for gRNA_sign in gRNA_sign_list:
        gRNA_count += 1
        gRNA = gRNAs[gRNA_sign]
        # generate sub cassette

        sub_cassette_sign = []
        sub_cassette = ""

        if starting_fragment:
            sub_cassette_sign.append("MtU6_3'20bp")
            sub_cassette = sub_cassette + product_start
            starting_fragment = False
        else:
            (sign, sequence) = draw_random_element(promoters)
            sub_cassette_sign.append(sign)
            sub_cassette = sub_cassette + sequence

        sub_cassette_sign.append(gRNA_sign)
        sub_cassette = sub_cassette + gRNA

        if gRNA_count == gRNA_num:
            # get the ending part
            sub_cassette_sign.append("SF_5'20bp")
            sub_cassette = sub_cassette + product_end
        else:
            # keep extending randomly
            (sign, sequence) = draw_random_element(scaffolds)
            sub_cassette_sign.append(sign)
            sub_cassette = sub_cassette + sequence

            (sign, sequence) = draw_random_element(linkers)
            sub_cassette_sign.append(sign)
            sub_cassette = sub_cassette + sequence

        # put the cassette to synthesis fragment and hypothetical product
        current_synthesis_fragment_sign.extend(sub_cassette_sign)
        current_synthesis_fragment = current_synthesis_fragment + sub_cassette
        hypothetical_product_sign.extend(sub_cassette_sign)
        hypothetical_product = hypothetical_product + sub_cassette

        # finishing or extending fragment
        if len(current_synthesis_fragment_sign) > fragment_element_length:  # finishing

            # store the finished fragment
            synthesis_fragments_sign.append(current_synthesis_fragment_sign)
            synthesis_fragments.append(current_synthesis_fragment)

            # re init current_synthesis_fragment with previous ending linker
            current_synthesis_fragment_sign = [sign]
            current_synthesis_fragment = sequence

    # check the basic sanity for each synthesis_fragment
    quality_control_list = [condition_check(
        synthesis_fragment) for synthesis_fragment in synthesis_fragments]

    if all(quality_control_list):
        product_count += 1
        for i in range(len(synthesis_fragments)):
            fragment_count += 1
            print(
                f'progress: {product_count}/{product_number}; {simulation_count} simulated')
            product_name = f'{run_name}_{product_count}/{simulation_count}_p{fragment_count}_{"-".join(synthesis_fragments_sign[i])}'
            # print(f'>{product_name}')
            # print(synthesis_fragments[i])
            product_dict[product_name] = synthesis_fragments[i]
            seq = synthesis_fragments[i]
            length = len(seq)
            (mfe_struct, mfe) = get_structure(seq)

            serial_number_list.append(
                f"{run_name}_{product_count}/{simulation_count}")
            part_list.append(f"p{fragment_count}")
            content_sign_list.append(
                f'{"-".join(synthesis_fragments_sign[i])}')
            seq_list.append(synthesis_fragments[i])
            length_list.append(length)
            mfe_list.append(mfe)
            mfe_struct_list.append(mfe_struct)

        if product_count == product_number:
            simulating = False
    if halt == simulation_count:
        print("halt")
        simulating = False

generated_table = pd.DataFrame.from_dict(
    {
        "serial_number": serial_number_list,
        "part": part_list,
        "content_sign": content_sign_list,
        "seq": seq_list,
        "length": length_list,
        "mfe": mfe_list,
        "mfe_struct": mfe_struct_list,
        # "define_structure":define_structure_list
    }
)
# TODO remove the identical set
generated_table = generated_table.join(generated_table.groupby(
    'serial_number')['mfe'].mean(), on='serial_number', rsuffix='_mean')
generated_table = generated_table.join(generated_table.groupby(
    'serial_number')['mfe'].std(), on='serial_number', rsuffix='_std')

generated_table = generated_table.sort_values(by='mfe_std', ascending=True)
generated_table = generated_table.sort_values(by='mfe_mean', ascending=False)

generated_table.to_csv(f"{run_name}_analysis_generated_seq.csv")

sns_plot = sns.scatterplot(data=generated_table, x='mfe_std', y='mfe')
sns_plot.figure.savefig(f"{run_name}_analysis_generated_seq.png")

with open(f"{run_name}_analysis_generated_seq.fa", 'w') as f:
    for key, value in product_dict.items():
        f.write(f'>{key}\n')
        f.write(f'{value}\n')
