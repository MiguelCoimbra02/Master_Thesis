// prevents dialog closing immediately when page navigates
vex.defaultOptions.closeAllOnPopState = false;

var netviz = {
    nodes: undefined,
    edges: undefined,
    network: undefined,
    isFrozen: false,
    newNodes: undefined,
    newEdges: undefined
};


// var network = null;
var node_search_data = null;
var node_search_data_dict = null;
var select = null;
var selectedRanks = new Set([0, 4]);
var rangeSliderValue = 0
var rangeSliderValue_viz_nodes = 100

let validNodes = []

//format.extend (String.prototype, {});

$(window).resize(function() {
    scale();
});

function enableSpinner() {
    // disable button
    $('#searchButton').prop("disabled", true);
    // add spinner to button
    $('#searchButton').html(`<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> searching...`);
}

function disableSpinner() {
    // enable button
    $('#searchButton').prop("disabled", false);
    // add back text only
    $('#searchButton').html('search');
}

// ajax request once the html is loaded to get node info
$(document).ready(function() {
    // Fetch node data from the server

    $('#dropdownMenu a').click(function(event){
        event.preventDefault(); // Prevent default page scroll when pressing <a> with tag href='#

        if ($(this).attr('href') == '#nodes') {
            export_nodes();
        }
        else if ($(this).attr('href') == '#edges') {
            export_edges();
        } else if ($(this).attr('href') == '#png') {
            export_network_png();
        } else if ($(this).attr('href') == '#report') {
            generate_report(); 
            // Handle user guide file download ->
        } else if ($(this).attr('href') == '#connecTF') {
            downloadFile('QS_output_final.csv');
        } else if ($(this).attr('href') == '#cis') {
            downloadFile('ALL_TF_Matches_Output.txt'); 
        } else if ($(this).attr('href') == '#dap') {
            downloadFile('dap_seq.txt');
        } else if ($(this).attr('href') == '#all_nodes') {
            downloadFile('nodes.csv');
        } else if ($(this).attr('href') == '#all_edges') {
            downloadFile('edges.csv');
        } else if ($(this).attr('href') == '#co_exp') {
            downloadFile('network.csv');
        }
    });

    function downloadFile(filename) {
        console.log("Requesting download for:", filename);
        fetch('/downloads/' + filename)
            .then(response => {
                if (response.ok) {
                    return response.blob(); // Get the file content
                }
                throw new Error('File not found');
            })
            .then(blob => {
                const url = URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = filename; // Set the file name for download
                document.body.appendChild(a);
                a.click(); // Trigger the download
                document.body.removeChild(a); // Clean up
            })
            .catch(error => {
                console.error('Error downloading file:', error);
            });
    }
    


    // Event listener for the search button
    $('#searchButton').click(function(event) {
        event.preventDefault();  // Prevent the default form submission and page reload
        console.log('ranks', selectedRanks)
        const userInput = $('input[name="query"]').val().trim();
        const queryNodes = userInput.split(/[,\s]+/);

        
        // If nothing, '', was searched show an error message
        if (queryNodes.length === 0) {
            alert('No gene was provided to do the search. Results showed are based on the ranks and Limited.');
            return;
        }
        enableSpinner();  // Show loading spinner
        
        // Make AJAX POST request
        $.ajax({
          url: "/search",
          dataType: 'json',
          type: "POST",
          contentType: 'application/json; charset=utf-8',
          processData: false,
          data: JSON.stringify({'nodes': Array.from(queryNodes),
                                'ranks': Array.from(selectedRanks),  // Convert Set to Array before sending
                                'rangeSliderValue': parseFloat(rangeSliderValue),
                                'rangeSliderValue_viz_nodes': parseFloat(rangeSliderValue_viz_nodes)

          }),
          success: function( data, textStatus, jQxhr ){
            disableSpinner();  // Hide loading spinner

            // Check if there are not found genes
            if (data.not_found && data.not_found.length > 0) {
                alert(`The following genes were not found in the database: ${data.not_found.join(', ')}`);
            } else if (data.not_found && data.not_found.length === 0 && (!data.found || data.found.length === 0)) {
                // If not_found is empty and found is also empty
                alert("None of the entered genes were found in the database.");
            }
            // Check if there are found genes and attribute them to the global variable
            
            if (data.max_result) {
                alert("The Limit of 100 nodes was achieved. Output was truncated.");
            }
            if (data.found) {
                validNodes = data.found;
            }
            
            if (data.graph) {
                netviz.isFrozen = false;
                drawNetwork(data.graph);
                console.log('draw_network', data.graph)
            }   
          },

          error: function(jqXhr) {
            if (jqXhr.responseJSON && jqXhr.responseJSON.error) {
                alert(jqXhr.responseJSON.error);
            } else {
                alert('Server error while loading the network.');
            }
          }
    });
    scale();
    initContextMenus();   
});


function drawNetwork(graphData){
    console.log("Graph Data Received:", graphData);
    console.log('edge data before viz', graphData.network.edges);

    netviz.nodes = new vis.DataSet(graphData.network.nodes);
    netviz.edges = new vis.DataSet(graphData.network.edges);
    // create a network
    var container = document.getElementById('networkView');

        // provide the data in the vis format
    var data = {
        nodes: netviz.nodes,
        edges: netviz.edges
    };
    var totalNodes = netviz.nodes.length;
    var totalEdges = netviz.edges.length;
    updateNodeCount(totalNodes);
    updateEdgeCount(totalEdges);

    updateNodeCount_TF_TR(netviz.nodes)
    updateDirect_Edges_count(netviz.edges)

    console.log("netviz.nodes:", netviz.nodes);
    console.log('netviz data', data);
    console.log('size', getNodeSize(totalNodes));
    console.log('total  nodes', totalNodes);
    console.log('total  edges', totalEdges);

    //https://visjs.github.io/vis-network/docs/network/#options

    
    var options = {height: '90%',
        interaction: {
            navigationButtons: true,
            keyboard: true,
            hover: true,
            multiselect:true,
                        
                    },
                    edges: {
                        //arrows: 'to',
                        smooth: {
                            enabled: true,
                            type: 'dynamic',
                            forceDirection: 'none'
                        },
                        font: {
                            size: 12,
                            face: 'sans',
                            align: 'top',
                            color: '#332f2f',
                        },
                        chosen: {
                            label: hover_edge_label
                        },
                        color: {color: '#6d468f', hover: '#513cd6', highlight: '#513cd6'},
                        hoverWidth: 2.2,
                        width: 2
                        },
                    nodes: {
                        shape: 'hexagon',
                        color: '#c7bba9',
                        size: getNodeSize(totalNodes),
                        font: {
                            multi: 'html'
                        },
                        chosen: {
                            node: hover_node,
                            label: hover_node_label
                        }
                    },  
                    physics: {
                        enabled: true,
                        solver: 'barnesHut',
                        barnesHut: {
                            gravitationalConstant: -25000,
                            centralGravity: 0.01,
                            springLength: 150,
                            springConstant: 0.08,
                            damping:  netviz.nodes.length > 100 ? 0.5 : 0.2,
                        },
                        repulsion: {
                            centralGravity: 0,
                            springLength: 200,
                            springConstant: 0.05,
                            nodeDistance: 150,
                            damping: 0.1
                        },
                        stabilization: {
                                enabled: true,
                                iterations: 1000,
                                fit: true,
                                updateInterval: 50,
                            },
                    },

                    configure: {
                        enabled: false
                    },
                    layout :{
                        improvedLayout: totalNodes <= 200
                    }   
    };
    postprocess_edges(data.edges);
    postprocess_nodes(data.nodes);
    console.log('edge items',data.edges)
    console.log('node items',data.nodes)

    netviz.network = new vis.Network(container, data, options);
    netviz.network.on('dragStart', onDragStart);
    netviz.network.on('dragEnd', onDragEnd);
    netviz.network.setOptions({interaction:{tooltipDelay:360000}}); //disable info table to appear when hovered
    netviz.network.on("stabilizationIterationsDone", function (params) {
        disableSpinner();
    });

    }
}); 

function hover_edge_label(values, id, selected, hovering) {
  values.mod = 'normal';
}

function hover_node_label(values, id, selected, hovering) {
  values.mod = 'normal';
}

function hover_node(values, id, selected, hovering) {
  values.borderWidth = 2.2;
  values.borderColor = '#513cd6'
  // values.color = 'blue'
}




function postprocess_edge(item) {

    let maxlen = 100;
    let header = '<div style="display: flex; justify-content: center; width: 100%;">\
                  <div style="overflow-x: auto; min-width: 700px; max-width: 1100px; width: 100%;">\
                  <table class="table table-striped table-bordered tooltip_table" style="width: 100%;">\
                  <tbody>';

    let footer = '</tbody></table>';
    let data = [];
    if (item.directed == '1') {
        data = [
            ['Edge ID', item.id],
            ['TF Rank', item.tf_rank],
            ['ConnecTF Analysis', item.ConnecTF],
            ['Cis-Elements analysis', item.cis_elements],
            ['Dap Seq Analysis', item.dap_seq],
            ['ConnecTF Method', item.edge_type],
            ['ConnecTF Predicted Product', item.gene_name],
            ['Cis-Elements Confidence value', item.cis_value],
            ['directed', item.directed]
        ];
    } else {
        data = [
            ['id', item.id],
            ['TF Rank', item.tf_rank],
            ['irp_score', item.irp_score],
            ['directed', item.directed]
        ];
    }

    // Generate the table rows
    let table = '';
    data.forEach(function (item, index) {
        if (item[1] != null) {
            let row = '<tr>\
                            <td><strong>' + item[0] + '</strong></td>\
                            <td class="text-wrap">' + item[1] + '</td>\
                       </tr>';
            table += row;
        }
    });
    table = header + table + footer;
    item.title = htmlTitle(table);
    return item;


}


function postprocess_edges(edges) {
    edges.forEach((item, i) => {
        edges[i] = postprocess_edge(item);
    });
}

function postprocess_node(node) {
    
    let maxlen = 100;
    

    let header = '<div style="display: flex; justify-content: center; width: 100%;">\
                <div style="overflow-x: auto; min-width: 700px; max-width: 1100px; width: 3000px%;">\
                <table class="table table-striped table-bordered tooltip_table" style="width: 100%;">\
                <tbody>';

    let footer = '</tbody>\
                  </table>\
                  </div>';

    let data = [['Node Name', node.name],
                ['homolog of Arabidopsis Concise', node.arabidopsis_gene],
                ['isTF', node.isTF == '1' ? 'True' : 'False'],
                ['isTR', node.isTR == '1' ? 'True' : 'False'],
                ['Gene Annotation', node.gene_annotation]];

    let table = '';
    data.forEach(function (pair) {
        if (pair[1] != null) {
            let row = '<tr>\
                            <td><strong>' + pair[0] + '</strong></td>\
                            <td class="text-wrap">' + pair[1] + '</td>\
                       </tr>';
            table += row;
        }
    });
    table = header + table + footer;
    node.title = htmlTitle(table);
    return node;
}

function postprocess_nodes(nodes) {
    nodes.forEach((item, i) => {
        nodes[i] = postprocess_node(item);
    });
}

function htmlTitle(html) {
    const container = document.createElement("div");
    container.classList.add('node_tooltip');
    container.style.display = "flex";
    container.style.justifyContent = "center";  // Center horizontally
    container.style.alignItems = "center";  // Align vertically
    container.style.width = "100%";  // Make sure the tooltip expands properly
    container.style.maxWidth = "900px";  // Limit maximum width
    container.style.whiteSpace = "normal";  
    container.style.padding = "10px";  
    container.style.border = "1px solid #ccc";  
    container.style.background = "#fff";  
    container.innerHTML = html;
    return container;
}
  




function scale() {
    $('#networkView').height(verge.viewportH());
    $('#networkView').width($('#networkViewContainer').width());
}

function freezeNodes(state){
    // Allows, or not, to move the nodes
    netviz.network.setOptions( { physics: !state } );
}

function onDragStart(obj) {
    if (obj.hasOwnProperty('nodes') && obj.nodes.length==1) {
        var nid = obj.nodes[0];
        netviz.nodes.update({id: nid, fixed: false});
    }

}

function onDragEnd(obj) {
    if (netviz.isFrozen==false)
        return
    var nid = obj.nodes;
    if (obj.hasOwnProperty('nodes') && obj.nodes.length==1) {
        var nid = obj.nodes[0];
        netviz.nodes.update({id: nid, fixed: true});
    }
}

function formatNodeInfoVex(nid) {
    return netviz.nodes.get(nid).title;
}

function formatEdgeInfoVex(nid) {
    return netviz.edges.get(nid).title;
}

function edge_present(edges, newEdge) {
    var is_present = false;
    var BreakException = {};

    try {
        edges.forEach((oldEdge, i) => {
            if (newEdge.id == oldEdge.id) {
                    is_present = true;
                    throw BreakException; // break is not available in forEach
                }
        })
    } catch (e) {
        if (e !== BreakException) throw e;
    }
    return is_present;
}

function expandNode(nid) {
    $.ajax({
      url: "/expand",
      async: false,
      dataType: 'json',
      type: "POST",
      contentType: 'application/json; charset=utf-8',
      processData: false,
      data: JSON.stringify({'node_name': nid, 
                            //'all_nodes': netviz.nodes.getIds(),
                            'ranks': Array.from(selectedRanks),  // Convert Set to Array before sending
                            'rangeSliderValue': parseFloat(rangeSliderValue)
                        }),
      success: function( data, textStatus, jQxhr ){
        console.log('expandnode data sent to app',data)
        if (data.error) {
            vex.dialog.alert('Server error when expanding the node.')
        }
        else {
            let nodes = data.graph.network.nodes;
            let edges = data.graph.network.edges;

            let newCounter = 0;

            if (data.max_result) {
                alert("The Limit of 100 nodes was achieved. Output was truncated.");
            }

            nodes.forEach((item) => {
                if (!netviz.nodes.get(item.id)) {
                    postprocess_node(item);
                    netviz.nodes.add(item);
                    newCounter += 1;
                }
            });
            
            edges.forEach((newEdge) => {
                if (!edge_present(netviz.edges, newEdge)) {
                    postprocess_edge(newEdge);
                    netviz.edges.add(newEdge);
                    newCounter += 1;
                }
            });

            if (newCounter==0) {
                vex.dialog.alert('No nodes or edges can be added.');
            }

            updateNodeCount(netviz.nodes.length);
            updateEdgeCount(netviz.edges.length);
            updateNodeCount_TF_TR(netviz.nodes);
            updateDirect_Edges_count(netviz.edges);
        }
        
    },
    error: function( jqXhr, textStatus, errorThrown ){
        alert('Server error while loading node data.');
    }
  });

}

function initContextMenus() {
    var canvasMenu = {
        "stop": {name: "Stop simulation"},
        "start" : {name: "Start simulation"}
    };
    var canvasMenu = {
        "freeze": {name: "Freeze positions"},
        // "release" : {name: "Start simulation"}
    };
    var nodeMenuFix = {
        "delete": {name: "Delete"},
        "expand": {name: "Expand"},
        "fix": {name: "Fix position"},
        "info": {name: "Info"}
    };
    var nodeMenuRelease = {
        "delete": {name: "Delete"},
        "expand": {name: "Expand"},
        "release": {name: "Release position"},
        "info": {name: "Info"}
    };
    var edgeMenu = {
        "delete": {name: "Delete"},
        "info": {name: "Info"}
    };

    $.contextMenu({
        selector: 'canvas',
        build: function($trigger, e) {
            // this callback is executed every time the menu is to be shown
            // its results are destroyed every time the menu is hidden
            // e is the original contextmenu event, containing e.pageX and e.pageY (amongst other data)

            var hoveredEdge = undefined;
            var hoveredNode = undefined;
            if (!$.isEmptyObject(netviz.network.selectionHandler.hoverObj.nodes)) {
                hoveredNode = netviz.network.selectionHandler.hoverObj.nodes[Object.keys(netviz.network.selectionHandler.hoverObj.nodes)[0]];
            }
            else {
                hoveredNode = undefined;
            }
            if (!$.isEmptyObject(netviz.network.selectionHandler.hoverObj.edges)) {
                hoveredEdge = netviz.network.selectionHandler.hoverObj.edges[Object.keys(netviz.network.selectionHandler.hoverObj.edges)[0]];
            }
            else {
                hoveredEdge = undefined;
            }

            // ignore auto-highlighted edge(s) on node hover
            if (hoveredNode != undefined && hoveredEdge != undefined)
                hoveredEdge = undefined;

            if (hoveredNode != undefined && hoveredEdge == undefined) {
                return {
                    callback: function(key, options) {
                        if (key == "delete") {
                            netviz.nodes.remove(hoveredNode);
                        }
                        else if (key == "expand") {
                            console.log(hoveredNode.id)
                            expandNode(hoveredNode.id);
                            //vex.dialog.alert("Not yet implemented.");
                        }
                        else if (key == "fix") {
                            netviz.nodes.update({id: hoveredNode.id, fixed: true});
                        }
                        else if (key == "release") {
                            netviz.nodes.update({id: hoveredNode.id, fixed: false});
                        }
                        else if (key == "info") {
                            vex.dialog.alert({unsafeMessage: formatNodeInfoVex(hoveredNode.id)});
                        }
                    },
                    items: netviz.nodes.get(hoveredNode.id).fixed ? nodeMenuRelease : nodeMenuFix
                };
            }
            else if (hoveredNode == undefined && hoveredEdge != undefined) {
                return {
                    callback: function(key, options) {
                        if (key == "delete") {
                            netviz.edges.remove(hoveredEdge);
                        }
                        else if (key == "info") {
                            vex.dialog.alert({unsafeMessage: formatEdgeInfoVex(hoveredEdge.id)});
                        }
                    },
                    items: edgeMenu
                };
            }
            else {
                if (netviz.isFrozen) {
                    canvasMenu.freeze.name = "Release positions";
                    return {
                        callback: function(key, options) {
                            if (key == "freeze") {
                                netviz.isFrozen = false;
                                freezeNodes(netviz.isFrozen);
                            }
                        },
                        items: canvasMenu
                    };
                }
                else {
                    canvasMenu.freeze.name = "Freeze positions";
                    return {
                        callback: function(key, options) {
                            if (key == "freeze") {
                                netviz.isFrozen = true;
                                freezeNodes(netviz.isFrozen);
                            }
                        },
                        items: canvasMenu
                    };
                }
            }
        }
    });

}


function format_cell(s){
    s = s.toString();
    s = s.trim();
    s = s.replace('\n', '');
    if (s[0]!='"' && s.slice(-1)!='"' && s.search(',')!=-1){
        s = '"' + s + '"';
    }
    return s;
}


// Toggle dropdown visibility
function toggleDropdown() {
    document.getElementById("dropdownContent").classList.toggle("show");
}

function toggleRank(rank) {
    let regulatoryCheckbox = document.querySelector('input[value="4"]'); // Select the Regulatory edges checkbox
    let selectedCheckbox = document.querySelector(`input[value="${rank}"]`);

    if (selectedCheckbox.checked) {
        // If the selected rank is 1, 2, or 3, visually select Regulatory edges
        if ([1, 2, 3].includes(rank)) {
            regulatoryCheckbox.checked = true; // Visually select Regulatory edges
            selectedRanks.add(4);
        }
        selectedRanks.add(rank); 
    } else {
        selectedRanks.delete(rank); 
    }

    console.log(Array.from(selectedRanks)); // Log selected ranks for debugging
}

//range slider
function updateRangeValue(value) {
    // Ensure that the value has 3 decimal places
    document.getElementById('manualInput').value = parseFloat(value).toFixed(3); // Update input with slider value
    // Sync the value with the slider
    document.getElementById('irpScoreRange').value = value; 
}

// Sync the slider with the manual input field when the user types a value
function syncSliderWithInput(value) {
    // Ensure the value is within range
    if (value < 0 || value > 1) {
        alert('Please enter a value between 0 and 1.');
        
    }

    // Store the cursor position
    let inputField = document.getElementById('manualInput');
    let cursorPosition = inputField.selectionStart;

    // Update slider value
    document.getElementById('irpScoreRange').value = value;

    // Only update if needed to avoid resetting cursor position
    if (inputField.value !== value) {
        inputField.value = parseFloat(value).toFixed(3);
    }
}

// Function to filter elements based on the selected IRP score
function filterByIrpScore(minIrpScore = 0) {
    console.log("Filtering edges with IRP Score >= " + minIrpScore);
    rangeSliderValue = minIrpScore
}

// slide ranger for limit_viz_node
function syncSliderWithInput_viz_nodes(value) {
    // Ensure the value is within range
    if (value < 0 || value > 300) {
        alert('Please enter a value between 0 and 300. The default is 100.');
        
    }

    // Store the cursor position
    let inputField = document.getElementById('manualInput_viz_nodes');
    let cursorPosition = inputField.selectionStart;

    // Update slider value
    document.getElementById('irpScoreRange_viz_nodes').value = value;

    // Only update if needed to avoid resetting cursor position
    if (inputField.value !== value) {
        inputField.value = parseFloat(value);
    }
}

// Function to filter elements based on the selected IRP score
function limit_viz_nodes(num_viz_nodes = 0) {
    console.log("Number of nodes showed in the vizualisation, after seacrhing, is " + num_viz_nodes);
    rangeSliderValue_viz_nodes = num_viz_nodes
}



function export_nodes() {
    console.log("Export nodes function called.");

    $.ajax({
        url: "/export_nodes",
        method: "GET",
        xhrFields: {
            responseType: 'blob'
        },
        success: function (data, status, xhr) {
            var filename = "nodes.csv";
            var disposition = xhr.getResponseHeader('Content-Disposition');
            if (disposition && disposition.indexOf('filename=') !== -1) {
                filename = disposition.split('filename=')[1].replace(/['"]/g, '');
            }

            var url = window.URL.createObjectURL(data);
            var a = document.createElement('a');
            a.href = url;
            a.download = filename;
            document.body.appendChild(a);
            a.click();
            a.remove();
        },
        error: function (xhr, status, error) {
            console.error("Error exporting nodes CSV:", error);
            vex.dialog.alert('Error exporting nodes CSV. Check the console for details.');
        }
    });
}


function export_edges() {
    console.log("Export edges function called.");

    $.ajax({
        url: "/export_edges",
        method: "GET",
        xhrFields: {
            responseType: 'blob'  // important: tells AJAX this is a binary file
        },
        success: function (data, status, xhr) {
            // Try to extract the filename from the Content-Disposition header
            var filename = "edges.csv";
            var disposition = xhr.getResponseHeader('Content-Disposition');
            if (disposition && disposition.indexOf('filename=') !== -1) {
                filename = disposition.split('filename=')[1].replace(/['"]/g, '');
            }

            var url = window.URL.createObjectURL(data);
            var a = document.createElement('a');
            a.href = url;
            a.download = filename;
            document.body.appendChild(a);
            a.click();
            a.remove();
        },
        error: function (xhr, status, error) {
            console.error("Error exporting CSV:", error);
            vex.dialog.alert('Error exporting CSV. Check the console for details.');
        }
    });
}


function export_network_png() {
    
    if (!netviz || !netviz.network) {
        vex.dialog.alert('No network to export! You need to do a search first.');
        return;
    }
    netviz.network.fit();
    // Give some time for the zoom effect to complete before capturing the canvas
    setTimeout(function() {
        // Get the canvas element inside the div with class 'vis_network'
        
        var canvas = document.querySelector('.vis-network canvas');
        
        if (canvas) {
            // Convert the canvas to a data URL (base64 encoded PNG)
            var dataURL = canvas.toDataURL("image/png");
            // Create a temporary link element to download the PNG
            var link = document.createElement('a');
            link.href = dataURL;
            link.download = 'network.png';  // Filename of the exported PNG
            document.body.appendChild(link); // Append link to the body
            link.click(); // Trigger download
            document.body.removeChild(link); // Remove the link after download
        } else {
            vex.dialog.alert('No network canvas found!');
            return;
        }
    }, 500);
}

function generate_report() {
    if (netviz.edges === undefined || netviz.nodes === undefined) {
        vex.dialog.alert('No data to create a report yet! You need to do a search first.');
        return;
    }
    
    console.log('Current netviz.edges:', netviz.edges.get());
    console.log('netviz.nodes', netviz.nodes.get());
    
    netviz.network.fit();
    
    setTimeout(function() {
        // Create new dictionaries to store only edge and node data
        const edge_data_dict = {};
        const node_data_dict = {};

        // Loop through the nodes using the get() method
        netviz.nodes.get().forEach(node => {
            // Check if the node object contains 'isTF' and 'id'
            if (node && typeof node === 'object' && node.hasOwnProperty('isTF') && node.hasOwnProperty('id')) {
                // Add the node's id as the key and the node itself as the value to the dictionary
                node_data_dict[node.id] = node;
            }
        });

        // Loop through the edges using the get() method
        netviz.edges.get().forEach(edge => {
            // Check if the edge object contains 'from' and 'to'
            if (edge && typeof edge === 'object' && edge.hasOwnProperty('from') && edge.hasOwnProperty('to')) {
                // Add the edge's id as the key and the edge itself as the value to the dictionary
                edge_data_dict[edge.id] = edge;
            }
        });

        // Capture the network image
        var canvas = document.querySelector('.vis-network canvas');
        var dataURL = canvas.toDataURL("image/png");

        // Send the data to the server
        $.ajax({
            url: '/report',
            type: 'POST',
            contentType: 'application/json; charset=utf-8',
            dataType: '',
            data: JSON.stringify({
                'quried_nodes': Array.from(validNodes),
                'tf_ranks': Array.from(selectedRanks),  // Convert Set to Array before sending
                'rangeSliderValue': parseFloat(rangeSliderValue),
                //'nodes': node_data_dict, iM getting this file from flask saves graph in search / expand
                //'edges': edge_data_dict, iM getting this file from flask saves graph in search / expand
                'network_image': dataURL
            }),
            xhrFields: {
                responseType: 'blob' // Set the response type to blob
            },
            success: function(data, textStatus, jQxhr) {
                disableSpinner();  // Hide loading spinner
                const url = window.URL.createObjectURL(data);
                const link = document.createElement('a');
                link.href = url;

                const baseFileName = 'gene_relation_report';
                const nodeNames = Array.from(validNodes).join('_'); 
                const pdfFileName = `${baseFileName}_${nodeNames}.pdf`;

                link.download = pdfFileName; // Set the dynamic file name
                document.body.appendChild(link);
                link.click();
                link.remove();
                window.URL.revokeObjectURL(url); // Clean up
            },
            error: function(jqXhr, textStatus, errorThrown) {
                disableSpinner();  // Hide loading spinner
                console.error('Error generating report:', errorThrown);
                vex.dialog.alert('Server error while generating the report.');
            }
        });
    }, 500);
}

function navToggleDropdown() {
    const dropdownMenu = document.getElementById("dropdownMenu");
    dropdownMenu.classList.toggle("show");
}

// Close the dropdown menu if the user clicks outside of it
window.onclick = function(event) {
    if (!event.target.matches('.dropdown-toggle')) {
        const dropdowns = document.getElementsByClassName("dropdown-menu");
        for (let i = 0; i < dropdowns.length; i++) {
            const openDropdown = dropdowns[i];
            if (openDropdown.classList.contains('show')) {
                openDropdown.classList.remove('show');
            }
        }
    }
};


// Function to calculate node size with 'var'
function getNodeSize(nodeCount) {
    var minNodes = 1;  // Minimum number of nodes
    var maxNodes = 75;  // Maximum number of nodes
    var minSize = 15;  // Minimum size of nodes
    var maxSize = 40;  // Maximum size of nodes

    if (nodeCount >= maxNodes) {
        return minSize;  // Fixed to return the minimum size for large graphs
    }
    
    if (nodeCount <= minNodes) {
        return maxSize;  // Return maximum size for smallest graph
    }
    
    // Correct scaling calculation
    var scalingFactor = (nodeCount - minNodes) / (maxNodes - minNodes);
    var size = maxSize - scalingFactor * (maxSize - minSize);
    return size;
}
function updateNodeCount(numberOfNodes) {
    document.getElementById('node-count').innerHTML = `<strong>${numberOfNodes}</strong>`;
}

function updateEdgeCount(numberOfEdges) {
    document.getElementById('edge-count').innerHTML = `<strong>${numberOfEdges}</strong>`;
}


function updateNodeCount_TF_TR(nodesDataSet) {
    let nodes = nodesDataSet.get(); // Convert DataSet to an array
    let numberOfTF = 0;
    let numberOfTR = 0;

    for (let i = 0; i < nodes.length; i++) {
        if (nodes[i].isTF === '1') {
            numberOfTF += 1;
        } else if (nodes[i].isTR === '1') {
            numberOfTR += 1;
        }
    }
    
    document.getElementById('tf-count').innerHTML = `<strong>${numberOfTF}</strong>`;
    document.getElementById('tr-count').innerHTML = `<strong>${numberOfTR}</strong>`;
}

function updateDirect_Edges_count (egdesDataSet) {
    let edges = egdesDataSet.get();
    let numberOfDirectEdges = 0;

    for (let i = 0; i < edges.length; i++) {
        if (edges[i] && edges[i].directed && edges[i].directed === '1') {
            numberOfDirectEdges += 1;
        }
    }
    // Update the content with the strong tag intact
    document.getElementById('direct-edge-count').innerHTML = `<strong>${numberOfDirectEdges}</strong>`;

}

/*
// filter page disappers if the click is not in the button or its components
window.onclick = function(event) {
    const dropdown = document.getElementById("dropdownContent");
    const button = document.querySelector('.dropbtn');

    // Check if the click was outside the dropdown and the button
    if (!event.target.matches('.dropbtn') && !dropdown.contains(event.target)) {
        if (dropdown.classList.contains('show')) {
            dropdown.classList.remove('show'); // Hide dropdown
        }
    }
}
*/


