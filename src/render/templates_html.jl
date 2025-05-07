# TODO make constant later

html_template = mt"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{{:protein_name}}} motifs</title>
    <link rel="stylesheet" href="styles.css">
</head>
<body>    
    <br><br>
    <div class="wrapper">

        <div id="nav" style="display: flex; justify-content: center;"></div>
        <br><br><br><br><br>
        <div class="container">
            <!-- Image Set 1 -->
            {{#:DF}}
            <div class="sliderGroup">
                <div class="imageTextContainer">
                    <div id={{:div_img_id}} class="imageContainer">
                        <img id="img{{:i}}" src="{{{:img_src}}}" alt="{{:img_alt}}">
                    </div>
                    <div id="{{:div_text_id}}" class="textContainer">
                        <p id="text{{:i}}_1" class="imageText">{{{:p_id1_default}}}</p>
                        <p id="text{{:i}}_2" class="imageText">{{{:p_id2_default}}}</p>
                        <p id="text{{:i}}_3" class="imageText">{{{:p_id3_default}}}</p>
                        <p id="text{{:i}}_4" class="imageText">{{{:p_id4_default}}}</p>
                        <p id="text{{:i}}_5" class="imageText">{{{:p_id5_default}}}</p>
                        <p id="text{{:i}}_6" class="imageText">{{{:p_id6_default}}}</p>
                    </div>
                </div>
                <div id={{:div_slide_id}} class="sliderContainer">
                    <!--  <span>Select Image:</span> -->
                    <input id="valR{{:i}}" type="range" min="0" max="{{:max_comb}}" value="0">
                    <span id="range{{:i}}">Image 1</span>
                </div>
            </div>
            {{/:DF}}
        </div>
    </div>
    <script src="scripts{{:j}}.js"></script>
    <!-- Modal Structure -->
    <div id="highlightModal" class="modal">
        <div class="modal-content">
            <span class="close" onclick="closeModal()">&times;</span>
            <div id="highlightedSequences"></div>
        </div>
    </div>

    <!-- Modal Structure for clustering -->
    <div id="highlightModal_cluster">
        <div id="highlightContent">
            <span class="close" onclick="closeModal_cluster()">&times;</span>
            <div class="modal-column">
                <img id="modalImage" src="" alt="Image">
            </div>
            <div class="modal-column">
                <p id="modalText">This is placeholder text. You can replace it with dynamic content as needed.</p>
            </div>
        </div>
    </div>

    <!-- Modal Structure for consensus_str -->
    <div id="highlightModal_text">
        <div id="highlightModal_text_content">
            <span class="close" onclick="closeModal_text()">&times;</span>
            <div class="modal-column">
                <p id="modalText1"> This is placeholder text. You can replace it with dynamic content as needed.</p>
                <button id="copyButton" onclick="copyText()">copy string</button>
            </div>
        </div>
    </div>

    <!-- Modal Structure for just image -->
    <div id="highlightModal_img">
        <div id="highlightModal_img_content">
            <span class="close" onclick="closeModal_img()">&times;</span>
            <div class="modal-column">
                <img id="modalImage1" src="" alt="Image">
            </div>
        </div>
    </div>
"""

# <!-- Modal Structure for clustering -->
#     <div id="highlightModal_cluster">
#         <div id="highlightContent">
#             <span class="close" onclick="closeModal_cluster()">&times;</span>
#             <div class="modal-column">
#                 <img id="modalImage" src="" alt="Image">
#             </div>
#             <div class="modal-column">
#                 <p id="modalText">This is placeholder text. You can replace it with dynamic content as needed.</p>
#             </div>
#         </div>
#     </div>



html_template_singleton = mt"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{{:protein_name}}} motifs</title>
    <link rel="stylesheet" href="styles.css">
</head>
<body>
    <br><br>
    <div class="wrapper">
        <div id="nav" style="display: flex; justify-content: center;"></div>
        <br><br><br><br><br>
        <div class="container">
            <!-- Image Set 1 -->
            {{#:DF}}
            <div class="sliderGroup">
                <div class="imageTextContainer">
                    <div id={{:div_img_id}} class="imageContainer">
                        <img id="img{{:i}}" src="{{{:img_src}}}" alt="{{:img_alt}}">
                        <br>
                        <span>{{:img_alt}}</span>
                    </div>
                    <div id="{{:div_text_id}}" class="textContainer">
                        <p id="text{{:i}}_1" class="imageText">{{{:p_id1_default}}}</p>
                        <p id="text{{:i}}_2" class="imageText">{{{:p_id2_default}}}</p>
                        <p id="text{{:i}}_3" class="imageText">{{{:p_id3_default}}}</p>
                        <p id="text{{:i}}_4" class="imageText">{{{:p_id4_default}}}</p>
                        <p id="text{{:i}}_5" class="imageText">{{{:p_id5_default}}}</p>
                    </div>
                </div>  
            </div>
            {{/:DF}}
        </div>
    </div>
    <script src="scripts{{:j}}.js"></script>

    <!-- Modal Structure -->
    <div id="highlightModal" class="modal">
        <div class="modal-content">
            <span class="close" onclick="closeModal()">&times;</span>
            <div id="highlightedSequences"></div>
        </div>
    </div>
    <!-- Modal Structure for just image -->
    <div id="highlightModal_img">
        <div id="highlightModal_img_content">
            <span class="close" onclick="closeModal_img()">&times;</span>
            <div class="modal-column">
                <img id="modalImage1" src="" alt="Image">
            </div>
        </div>
    </div>
"""

html_hover_default = mt"""
    <div {{{:hover_on}}}>
        <div class="hover-meta-data">
        {{{:meta_str}}}
        </div>
    </div>"""

html_end = "</body></html>";
