template_css = mt"""
body {
    font-family: Arial, sans-serif;
    background-color: #ffffff;
    margin: 0;
    padding: 0;
    justify-content: center;
    align-items: center;
    height: 100vh;
}

span.putBar {
  border-top: 1px solid #5144FA;;
}

.wrapper {
    width: 65%; /* Center 60% of the page */
    margin: 0 auto; /* Ensure it's centered */
}

.current {
    font-weight: bold;
    pointer-events: none;
    text-decoration: none; /* Remove underline */
}

.container {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(315px, 1fr)); /* Smaller columns for tighter layout */
    gap: 10px; /* Reduced spacing between items */
    justify-items: center;
}

.sliderGroup {
    margin-bottom: 40px;
}

.imageTextContainer {
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 5px;
}

.imageContainer {
    margin-right: 1px;
    font-size: 11px;
    align-items: center;
    justify-content: center;
    text-align: center; /* Center the text below the image */
    margin-right: -25px; 
}

.imageContainer img {
    width: 70%; /* Adjust the percentage as needed */
    height: auto;
    border-radius: 5px;
    /* transition: opacity 0.01s ease-in-out; */
}

.textContainer {
    width: auto; /* Reduced width to fit the text */
    text-align: center; /* Center the text below the image */
}

/* Keep special styling for the current page link */
.textContainer a {
    text-decoration: none; /* Remove underline */
    color: #5144FA; /* Light blue color for links */
    margin: 0 5px; /* Even spacing between links */
}
    
.imageText {
    background-color: rgba(0, 0, 0, 0.05);
    padding: 6px;
    border-radius: 15px;
    font-size: 8px;
    white-space: nowrap; /* Ensure text doesn't wrap to fit the width */
    overflow: hidden; /* Hide any overflow */
    text-overflow: ellipsis; /* Add ellipsis if the text overflows */
}

.sliderContainer {
    display: flex;
    flex-direction: column;
    align-items: center;
    font-size: 10px;
}

.sliderContainer input[type="range"] {
    width: 25%; /* Shortened width of the sliders */
    margin-top: 10px;
    -webkit-appearance: none;
    appearance: none;
    height: 10px;
    background: #ddd;
    outline: none;
    opacity: 0.65;
    /* transition: opacity 0.01s ease; */
    border: 1.5px solid gray; /* Thin black border around the slider */
    border-radius: 6px;
}

.sliderContainer input[type="range"]:hover {
    opacity: 1;
}

.sliderContainer input[type="range"]::-webkit-slider-thumb {
    -webkit-appearance: none;
    appearance: none;
    width: 17px;
    height: 17px;
    border-radius: 60%;    
    background: silver;
    border: 1.5px solid black;    
    cursor: pointer;
    /* transition: background 0.01s ease; */
    z-index: 1;
}

.sliderContainer input[type="range"]::-moz-range-thumb {
    width: 20px;
    height: 20px;
    border-radius: 50%;
    background: silver;
    border: 2px solid black;
    cursor: pointer;
    /* transition: background 0.01s ease; */
}

 /* Modal Styles (for sequence substring highlight)*/

.column {
    display: flex;
    flex-direction: column;
    width: 48%; /* Adjust width as needed */
}

.modal {
    display: none; /* Hidden by default */
    position: fixed; /* Stay in place */
    z-index: 1000; /* Sit on top */
    left: 0;
    top: 0;
    width: 100%;
    height: 100%;
    overflow: auto; /* Enable scroll if needed */
    background-color: rgba(0,0,0,0.4); /* Black with opacity */
}

.modal-content {
    background-color: #fefefe;
    margin: 20% auto; /* 20% from the top and centered */
    padding: 25px;
    border: 2px solid #888;
    width: 80%; /* Could be more or less, depending on screen size */
    max-width: 985px; /* Limit max width */
    font-family: Arial, sans-serif;
    position: relative;
}

.close {
    color: #aaa;
    float: right;
    font-size: 28px;
    font-weight: bold;
    cursor: pointer;
}

.close:hover,
.close:focus {
    color: black;
    text-decoration: none;
    cursor: pointer;
}

.highlight {
    font-weight: bold;
    color: orange;
}

.highlight-comp {
    font-weight: bold;
    color: LightSteelBlue;
}

.sequence {
    font-family: monospace;
    white-space: pre-wrap; /* Preserve whitespace and wrap text */
    margin: 0;
    padding: 5px;
    text-align: left; /* Ensure text is aligned to the left */
}

.header {
    font-family: monospace;
    margin: 1px 0; /* Space above and below headers */
}

/* Modal cluster styles */
#highlightModal_cluster {
   display: none;
   position: fixed;
   z-index: 1;
   left: 0;
   top: 0;
   width: 100%;
   height: 100%;
   background-color: rgba(0, 0, 0, 0.5);
   overflow: auto;
}

#highlightModal_text {
   display: none;
   position: fixed;
   z-index: 1;
   left: 0;
   top: 0;
   width: 100%;
   height: 100%;
   background-color: rgba(0, 0, 0, 0.5);
   overflow: auto;

    /* Allow dynamic width adjustment */
   justify-content: center; 
   align-items: center; /* Center the modal */
}

#highlightModal_img {
   display: none;
   position: fixed;
   z-index: 1;
   left: 0;
   top: 0;
   width: 100%;
   height: 100%;
   background-color: rgba(0, 0, 0, 0.5);
   overflow: auto;
}


#highlightModal_text_content {
    background-color: white;
    margin: 15% auto;
    padding: 20px;
    border: 1px solid #888;
    text-align: center;
    max-width: 90%; /* Prevent it from exceeding 90% of the screen width */
    min-width: 200px; /* Ensure it doesn't shrink too much */
    word-wrap: break-word; /* Break long words if needed */
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2); /* Optional: add some style */
    border-radius: 8px; /* Optional: rounded corners */
    display: flex;
}


#highlightModal_img_content {
    background-color: white;
    margin: 15% auto;
    padding: 20px;
    border: 1px solid #888;
    text-align: center;
    max-width: 90%; /* Prevent it from exceeding 90% of the screen width */
    min-width: 200px; /* Ensure it doesn't shrink too much */
    word-wrap: break-word; /* Break long words if needed */
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2); /* Optional: add some style */
    border-radius: 8px; /* Optional: rounded corners */
    display: flex;
}


#highlightContent {
   background-color: white;
   margin: 15% auto;
   padding: 20px;
   border: 1px solid #888;
   width: 80%;
   text-align: center;
   display: flex;
   justify-content: space-around;
   align-items: center;
}

#modalText {
    font-size: 12px; /* Adjust this value as needed */
}

#modalText1 {
    font-size: 12px; /* Adjust this value as needed */
}

#copyButton {
    padding: 12px 24px;
    background-color: white;
    color: gray; /* gray text color */
    border: 2px solid lightgray; /* gray border */
    border-radius: 50px; /* Rounded, pill-like shape */
    cursor: pointer;
    font-size: 12px;
    font-weight: 500;
    margin-top: 15px;
    text-align: center;
    transition: all 0.3s ease; /* Smooth transition */
    box-shadow: 0 4px 10px rgba(0, 123, 255, 0.1); /* Subtle shadow */
}

#copyButton:hover {
    background-color: lightgray; /* Blue background on hover */
    color: black; /* White text on hover */
    border-color: lightgray; /* Darker border on hover */
    box-shadow: 0 6px 15px rgba(0, 123, 255, 0.2); /* Stronger shadow on hover */
}

#copyButton:focus {
    outline: none; /* Removes outline when focused */
    box-shadow: 0 0 5px rgba(0, 123, 255, 0.5); /* Focused glow effect */
}

.modal-column {
   flex: 1;
   padding: 10px;
}

.modal-column img {
   max-width: 100%;
}

.hover-window {
    position: fixed;
    left: 20px; /* Position from the left */
    top: 100px; /* Position from the top */
    width: 125px;
    padding: 10px;
    background-color: rgba(255, 255, 255, 0.9);
    border: 1px solid #ccc;
    border-radius: 5px;
    box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
    transition: transform 0.2s;
    font-size: 11px; /* Adjust this value as needed */
}

.hover-meta-data {    
    font-size: 8px; /* Adjust this value as needed */
}

.cl {
    list-style-type: none; /* Remove default list style */
    padding: 0; /* Remove padding */
    margin-bottom: 5px;
    display: flex; /* Use flex to align items */
    align-items: center; /* Center vertically */
}

.color-square {
    width: 15px; /* Width of the color square */
    height: 15px; /* Height of the color square */
    margin-right: 10px; /* Space between square and text */
    border-radius: 3px; /* Optional: rounded corners */
}

/* General navigation styling */
#nav {
    font-family: 'Arial', sans-serif; /* Use a clean and widely available font */
    font-size: 14px; /* Readable font size */
    background-color: #f9f9f9; /* Light background for contrast */
    border: 1px solid #ddd; /* Subtle border */
    border-radius: 12px; /* Rounded corners for a modern look */
    text-align: center; /* Center align the text */
    padding: 5px 10px; /* Add some padding for spacing */
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1); /* Subtle shadow for depth */
    max-width: fit-content; /* Limit the width of the navigation box */
    margin: 0 auto; /* Center the nav box within the page */
}

/* Styling for navigation links */
#nav a {
    text-decoration: none; /* Remove underline */
    color: #007BFF; /* Light blue color for links */
    margin: 0 5px; /* Even spacing between links */
}

/* Current page link styling */
#nav a.current {
    font-weight: bold; /* Highlight the current page */
    color: #000; /* Darker color for the current page */
}

/* The horizontal line */
.horizontal-line {
   border-top: 1px solid #ddd; /* Thin horizontal line */
   margin: 10px 0; /* Space above and below the line */
 }


"""
