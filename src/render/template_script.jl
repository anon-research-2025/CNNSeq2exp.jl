# TODO make this constant

script_template = mt"""
// Function to load images, labels, and texts from JSON

function createArray(num) {
    return Array.from({ length: num }, (_, i) => i + 1);
}

function updateNav() {
    const currentPage = 'index{{:j}}.html'; // Set this to the current page name
    const availablePages = createArray({{:upto}}); // List of available page numbers

    const navHtml = availablePages.map(num => {
        const page = `index${num}.html`;
        const label = num === 1 ? 'Singleton motifs' :
                      num === 2 ? 'Paired motifs' :
                      num === 3 ? 'Triplet motifs' :
                      num === 4 ? 'Quadruplet motifs' :
                      num === 5 ? 'Quintuplet motifs' :
                      'Sextuplet motifs';

        return `<a href="${page}" ${currentPage === page ? 'class="current"' : ''}>${label}</a>`;
    }).join(' &nbsp&nbsp | &nbsp&nbsp ');

    // Add the "readme" link at the end
    const readmePage = '../../../readme/readme.html';
    const readmeLink = `<a href="${readmePage}" ${currentPage === readmePage ? 'class="current"' : ''}>Readme</a>`;

    // Append the readme link to the existing navHtml
    document.getElementById('nav').innerHTML = navHtml + ' &nbsp&nbsp | &nbsp&nbsp ' + readmeLink;
}

window.onload = updateNav;

function loadData(modeCount, jsonData, callback) {
    const images = [];
    const labels = [];
    const texts = [];

    for (let modeIndex = 1; modeIndex <= modeCount; modeIndex++) {
        const modeImages = jsonData[`mode_${modeIndex}`].pwms;
        const modeLabels = jsonData[`mode_${modeIndex}`].labels;
        const modeTexts = jsonData[`mode_${modeIndex}`].texts;

        images.push(modeImages);
        labels.push(modeLabels);
        texts.push(modeTexts);
    }

    callback(images, labels, texts);
}

// Initialize the sliders with the loaded data
function initializeSliders(images, labels, texts) {
    images.forEach((modeImages, index) => {
        const modeLabels = labels[index];
        const modeTexts = texts[index];
        const sliderId = `valR${index + 1}`;
        const imgElementId = `img${index + 1}`;
        const rangeElementId = `range${index + 1}`;
        const textElementId1 = `text${index + 1}_1`;
        const textElementId2 = `text${index + 1}_2`;
        const textElementId3 = `text${index + 1}_3`;
        const textElementId4 = `text${index + 1}_4`;
        const textElementId5 = `text${index + 1}_5`;
        const textElementId6 = `text${index + 1}_6`;
        const sliderElement = document.getElementById(sliderId);
        
        function showVal(newVal) {
            const imgElement = document.getElementById(imgElementId);
            imgElement.style.opacity = 0;  // Start fading out

            setTimeout(function() {
                document.getElementById(rangeElementId).innerHTML = modeLabels[newVal];
                document.getElementById(textElementId1).innerHTML = modeTexts[newVal][0];
                document.getElementById(textElementId2).innerHTML = modeTexts[newVal][1];
                document.getElementById(textElementId3).innerHTML = modeTexts[newVal][2];
                document.getElementById(textElementId4).innerHTML = modeTexts[newVal][3];
                document.getElementById(textElementId5).innerHTML = modeTexts[newVal][4];
                document.getElementById(textElementId6).innerHTML = modeTexts[newVal][5];
                imgElement.src = modeImages[newVal];
                imgElement.style.opacity = 1;  // Fade back in
            }, 250);  // Wait for half the transition duration (0.5s / 2)
        }

        // Set initial values for each slider
        document.getElementById(rangeElementId).innerHTML = modeLabels[0];
        document.getElementById(textElementId1).innerHTML = modeTexts[0][0];
        document.getElementById(textElementId2).innerHTML = modeTexts[0][1];
        document.getElementById(textElementId3).innerHTML = modeTexts[0][2];
        document.getElementById(textElementId4).innerHTML = modeTexts[0][3];
        document.getElementById(textElementId5).innerHTML = modeTexts[0][4];
        document.getElementById(textElementId6).innerHTML = modeTexts[0][5];
        document.getElementById(imgElementId).src = modeImages[0];

        // Attach event listener to each slider
        document.getElementById(sliderId).addEventListener('input', function() {
            showVal(this.value);
        });

        if (modeImages.length == 1) {
            sliderElement.style.display = 'none';
            document.getElementById(rangeElementId).innerHTML = modeLabels[0];
        }

    });
}

// Variable to store the loaded JSON data
let jsonData = {};

// Function to load JSON data and store it in the jsonData variable
function loadJSONData(callback) {
    fetch('data{{:j}}.json')
        .then(response => response.json())
        .then(data => {
            jsonData = data;
            console.log('JSON data loaded:', jsonData);
            if (callback) callback();
        })
        .catch(error => {
            console.error('Error loading JSON data:', error);
        });
}

// Call the function to load the JSON data
loadJSONData(() => {
    // Ensure jsonData is loaded and then call loadData
    loadData({{:mode_counts}}, jsonData, initializeSliders);
});

// --------------------------------------- sequence hightlighting part  ---------------------------------------------

let scrollPosition = 0;

function openHtmlWindowImg(imageFile) {
    scrollPosition = window.pageYOffset || document.documentElement.scrollTop;
    // Set the image source in the modal
    const modalImage = document.getElementById('modalImage1');
    modalImage.src = imageFile;

    // Wait for the image to load to get its natural dimensions
    modalImage.onload = function () {
        // Get the image's natural width
        const imgWidth = modalImage.naturalWidth;

        // Dynamically set the modal width (optional max width for responsiveness)
        const modal = document.getElementById('highlightModal_img_content');
        modal.style.width = imgWidth > 800 ? '800px' : imgWidth + 'px'; // Cap at 800px for responsiveness
    };

    // Display the modal
    document.getElementById('highlightModal_img').style.display = "block";
}
   
function openHtmlWindowText(textContent) {
    // Store the current scroll position
    scrollPosition = window.pageYOffset || document.documentElement.scrollTop;

    // Set the text content
    const modalText = document.getElementById('modalText1');
    modalText.innerHTML = textContent;

    // Dynamically adjust the modal width based on text length
    const modalContent = document.getElementById('highlightModal_text_content');
    const contentLength = textContent.length;

    // Calculate width: 10px per character, with min and max limits
    const width = Math.min(Math.max(contentLength * 10, 200), 800); 
    modalContent.style.width = width + 'px';

    // Display the modal
    document.getElementById('highlightModal_text').style.display = "flex";
}

// Function to copy text from the modal
function copyText() {
    const modalTextElement = document.getElementById('modalText1');
    const originalText = modalTextElement.innerText;  // Save the original text
    const textToCopy = originalText;
    const textarea = document.createElement('textarea');
    textarea.value = textToCopy;
    document.body.appendChild(textarea);
    textarea.select();
    document.execCommand('copy');
    document.body.removeChild(textarea);

    // Display success message temporarily in the modal
    modalTextElement.innerHTML = 'string copied successfully!';

    // After 1 second, revert the content back to the original
    setTimeout(() => {
        modalTextElement.innerHTML = originalText;
    }, 1000);  // 1 second delay
}

function openHtmlWindow(imageFile, textContent) {
    // Store the current scroll position
    scrollPosition = window.pageYOffset || document.documentElement.scrollTop;

    // Set the image source in the modal
    document.getElementById('modalImage').src = imageFile;

    // Set the dynamic text content
    document.getElementById('modalText').innerHTML = textContent;

    // Display the modal
    document.getElementById('highlightModal_cluster').style.display = "block";
}

function loadFile(filePath) {
    return fetch(filePath).then(response => response.text());
}

function openHighlightPage(csvFile) {
    // Store the current scroll position
    scrollPosition = window.pageYOffset || document.documentElement.scrollTop;

    const fastaFiles = [
        {{#:sequence_file_paths}}
        '{{{.}}}'{{^.[end]}}, {{/.[end]}}
        {{/:sequence_file_paths}}
    ];

    // Load all FASTA files and the CSV file
    const fastaPromises = fastaFiles.map(file => loadFile(file));
    Promise.all([...fastaPromises, loadFile(csvFile)])
        .then(contents => {
            const fastaContents = contents.slice(0, fastaFiles.length);
            const csvContent = contents[contents.length - 1];

            // Combine sequences from all FASTA files
            const sequences = fastaContents.flatMap(fastaContent => parseFasta(fastaContent));
            const highlights = parseCsv(csvContent);

            const htmlContent = generateHtml(sequences, highlights);

            // Insert the generated HTML into the modal
            document.getElementById('highlightedSequences').innerHTML = htmlContent;

            // Show the modal
            document.getElementById('highlightModal').style.display = "block";
        });
}

function closeModal() {
    document.getElementById('highlightModal').style.display = "none";

    // Restore the scroll position
    window.scrollTo(0, scrollPosition);
}

function closeModal_text() {
    // Hide the modal
    document.getElementById('highlightModal_text').style.display = "none";

    // Restore scroll position
    window.scrollTo(0, scrollPosition);
}

function closeModal_cluster() {
    // Hide the modal
    document.getElementById('highlightModal_cluster').style.display = "none";

    // Restore scroll position
    window.scrollTo(0, scrollPosition);
}

function closeModal_img() {
    // Hide the modal
    document.getElementById('highlightModal_img').style.display = "none";

    // Restore scroll position
    window.scrollTo(0, scrollPosition);
}


// Close the modal when clicking outside of it
window.onclick = function(event) {
    const modal = document.getElementById('highlightModal');
    if (event.target === modal) {
        closeModal();
    }        
    const modal_cluster = document.getElementById('highlightModal_cluster');
    if (event.target === modal_cluster) {
        closeModal_cluster();
    }
    const modal_text = document.getElementById('highlightModal_text');
    if (event.target === modal_text) {
        closeModal_text();
    }
    const modal_img = document.getElementById('highlightModal_img');
    if (event.target === modal_img) {
        closeModal_img();
    }
}

// Close the modal when pressing the "Esc" key
window.onkeydown = function(event) {
    if (event.key === "Escape") {
        closeModal();
        closeModal_cluster();
        closeModal_text();
        closeModal_img();
    }
}

function parseFasta(fastaContent) {
    const lines = fastaContent.split('\n');
    const sequences = [];
    let currentHeader = '';
    let currentSequence = '';
    lines.forEach(line => {
        if (line.startsWith('>')) {
            if (currentHeader) {
                sequences.push({ header: currentHeader, sequence: currentSequence });
                currentSequence = '';
            }
            currentHeader = line;
        } else {
            currentSequence += line.trim();
        }
    });
    if (currentHeader) {
        sequences.push({ header: currentHeader, sequence: currentSequence });
    }
    return sequences;
}

function parseCsv(csvContent) {
    const lines = csvContent.trim().split('\n');
    const highlights = [];
    lines.slice(1).forEach(line => {
        const [seqIndex, startPosition, endPosition, iscomp] = line.split(',').map(Number);
        if (!isNaN(seqIndex) && !isNaN(startPosition) && !isNaN(endPosition) && !isNaN(iscomp)) {
            highlights.push({ seqIndex, startPosition, endPosition, iscomp });
        }
    });
    return highlights;
}

function generateHtml(sequences, highlights) {
    const uniqueIndices = [...new Set(highlights.map(h => h.seqIndex))];
    let htmlContent = '';

    uniqueIndices.forEach(index => {
        const sequence = sequences[index];
        const highlightsForSequence = highlights.filter(h => h.seqIndex === index);
        const highlightedSequence = highlightSequence(sequence.sequence, highlightsForSequence);

        htmlContent += `
            <div class="header">${sequence.header}</div>
            <div class="sequence">${highlightedSequence}</div>
        `;
    });

    return htmlContent;
}

function highlightSequence(sequence, highlights) {
    // Sort the highlights by startPosition
    highlights.sort((a, b) => a.startPosition - b.startPosition);

    // Merge overlapping intervals
    const mergedHighlights = [];
    highlights.forEach(({ startPosition, endPosition, iscomp }) => {
        if (
            mergedHighlights.length === 0 ||
            mergedHighlights[mergedHighlights.length - 1].endPosition < startPosition
        ) {
            mergedHighlights.push({ startPosition, endPosition, iscomp });
        } else {
            const last = mergedHighlights[mergedHighlights.length - 1];
            last.endPosition = Math.max(last.endPosition, endPosition);
            // If `iscomp` differs, we keep the last `iscomp` for simplicity
            if (last.iscomp !== iscomp) {
                last.iscomp = iscomp;
            }
        }
    });

    // Apply the highlights to the sequence
    let highlighted = '';
    let lastIndex = 0;
    
    mergedHighlights.forEach(({ startPosition, endPosition, iscomp }) => {
        if (startPosition > lastIndex) {
            highlighted += sequence.substring(lastIndex, startPosition);
        }

        // Apply a different color based on `iscomp`
        const colorClass = iscomp === 1 ? 'highlight-comp' : 'highlight';
        highlighted += `<span class="${colorClass}">${sequence.substring(startPosition, endPosition)}</span>`;
        lastIndex = endPosition;
    });

    if (lastIndex < sequence.length) {
        highlighted += sequence.substring(lastIndex);
    }


    return highlighted;
}

// hover window

    const hoverWindow = document.getElementById('hoverWindow');

    document.addEventListener('scroll', () => {
        const scrollY = window.scrollY || window.pageYOffset;
    });

"""
