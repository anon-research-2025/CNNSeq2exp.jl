const x_dir = [
    "singleton_motifs",
    "pairs_motifs",
    "triplet_motifs",
    "quadruplet_motifs",
    "quintuplet_motifs",
    "sextuplet_motifs"
]

const csv_header = [:seq_index, :start_position, :end_position, :is_reversed_complement]

const pwm_dpi = 65

const _ind2dna_str_ = Dict{Int, Char}(1 => 'A', 2 => 'C', 3 => 'G', 4 => 'T')
const _ind2dna_str_rna = Dict{Int, Char}(1 => 'A', 2 => 'C', 3 => 'G', 4 => 'U')
const _placeholder_char_ = 'n'

# html rendering
const pwms_str = "pwms"
const labels_str = "labels"
const texts_str = "texts"

# html tags
const tag_div_img_id = "div_img_id"
const tag_i = "i"
const tag_img_src = "img_src"
const tag_img_alt = "img_alt"
const tag_div_text_id = "div_text_id"
const tag_p_id1_default = "p_id1_default"
const tag_p_id2_default = "p_id2_default"
const tag_p_id3_default = "p_id3_default"
const tag_p_id4_default = "p_id4_default"
const tag_p_id5_default = "p_id5_default"
const tag_p_id6_default = "p_id6_default"
const tag_div_slide_id = "div_slide_id"
const tag_max_comb = "max_comb"
const DBSCAN_eps = 7.0

# boxplots
const row_width = 175
const col_width = 710
const text_scale_factor = 0.1
const boxplot_marker_size = 2