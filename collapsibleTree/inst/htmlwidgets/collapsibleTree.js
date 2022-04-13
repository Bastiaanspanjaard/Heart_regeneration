HTMLWidgets.widget({
  name: 'collapsibleTree',
  type: 'output',

  factory: function(el, width, height) {

  // PO define colors
  var fixed_colors = d3.scaleOrdinal().range(["#000023", "#03002C", "#070035", "#0A003E", "#0E0047", "#120050", "#150059", "#190063", "#1D006C", "#200075", "#24007E", "#280087", "#2B0090", "#2F0099", "#3300A3", "#3700AA", "#3C00AD", "#4200B1", "#4700B4", "#4D00B7", "#5200BA", "#5800BD", "#5D00C0", "#6300C3", "#6800C6", "#6E00C9", "#7300CD", "#7900D0", "#7E00D3", "#8400D6", "#8901D1", "#8F03C5", "#9406B8", "#9A08AC", "#A00A9F", "#A50D93", "#AB0F86", "#B0127A", "#B6146D", "#BB1761", "#C11955", "#C61C48", "#CC1E3C", "#D2202F", "#D72323", "#DA2E20", "#DD3D20", "#E04B20", "#E25920", "#E56720", "#E77620", "#EA8420", "#EC9220", "#EFA120", "#F2AF20", "#F4BD20", "#F7CC20", "#F9DA20", "#FCE820", "#FFF721"])


    var i = 0,
    duration = 750,
    root = {},
    options = {},
    treemap;

  // PO define pieNode size
     var pieNodeOut = 10,
	pieNodeInner= pieNodeOut/1e6,
	scSize = 4;
  // PO define pie and arc PO
  // TODO
     var arc = d3.arc().innerRadius(pieNodeInner),
         pie = d3.pie().sort(null);
    // Optionally enable zooming, and limit to 1/5x or 5x of the original viewport
    var zoom = d3.zoom()
    .scaleExtent([1/10, 10])
    .on('zoom', function () {
      if (options.zoomable) svg.attr('transform', d3.event.transform)
    })

    // create our tree object and bind it to the element
    // appends a 'group' element to 'svg'
    // moves the 'group' element to the top left margin
    var svg = d3.select(el).append('svg')
    .attr('width', width)
    .attr('height', height)
    .call(zoom)
    .append('g');

    // Define the div for the tooltip
    var tooltip = d3.select(el).append('div')
    .attr('class', 'tooltip')
    .style('opacity', 0);

    function update(source) {

      // PO If colors provided, use them
      var color;
      if(options.useColors){
	color = options.colors
	}else{
	color = fixed_colors
	}
      // Assigns the x and y position for the nodes
      var treeData = treemap(root);

      // Compute the new tree layout.
      var nodes = treeData.descendants(),
      links = treeData.descendants().slice(1);

      // Normalize for fixed-depth.
      nodes.forEach(function(d) {d.y = d.depth * options.linkLength});

      // ****************** Nodes section ***************************

      // Update the nodes...
      var node = svg.selectAll('g.node')
      .data(nodes, function(d) {return d.id || (d.id = ++i); });

      // Enter any new modes at the parent's previous position.
      var nodeEnter = node.enter().append('g')
      .attr('class', 'node')
       .attr('size', function(d) {
        return d.data.SizeOfNode;
      })
       .attr('isScar', function(d) {
        return (options.pieNode? d.data.isScar : false);
      })

      .attr('transform', function(d) {
        return 'translate(' + source.y0 + ',' + source.x0 + ')';
      })
      .on('click', click);

      // Add tooltips, if specified in options
      if (options.tooltip) {
        nodeEnter = nodeEnter
        .on('mouseover', mouseover)
        .on('mouseout', mouseout);
      }

      // PO append arcs
      if(options.pieNode){
       var arcGs = nodeEnter.selectAll("g.arc")
        .data(function(d) {
	//var proportions = (d.data.isScar?  d.data.pieNode : [1]) 
	var proportions = d.data.pieNode
	var prop_sum = proportions.reduce(function(a, b) { return a + b; }, 0);	
	if(prop_sum == 0){
		proportions = [1]
	}
	
	var angle = -JSON.parse(options.angle); //-Math.PI/2;  
	var pieProp = pie(proportions).map(function(m, j){
		m.id = j
		m.color = (prop_sum == 0 ? "#bbbbbb": color[j])
		m.isScar = JSON.parse(d.data.isScar);
		m.r = d.data.SizeOfNode;
		m.ct = (m.isScar? d.data.name : d.data.ct);
		m.startAngle = m.startAngle + angle;
		m.endAngle = m.endAngle + angle;
		return m;
		})
          return pieProp;
        });

        var arcEnter = arcGs.enter().append("g").attr("class", "arc");

        arcEnter.append("path")
         .attr("d", function(d) {
        size = 	(d.isScar ? d.r  : scSize)
         //arc.outerRadius(pieNodeOut);
         arc.outerRadius(size);
	var final_arc = arc(d);
         return final_arc;
         })
         .style("fill", function(d, i) { 
         return d.color });

      }

      // Add Circle for the nodes
      // PO conditional circles
      if(!options.pieNode){
       nodeEnter.append('circle')
       .attr('class', 'node')
       .attr('r', 1e-6)
       .style('fill', function(d) {
         //return d.data.fill || (d._children ? options.fill : color[1]);
         return d.data.fill || (d._children ? options.fill : '#fff');
       })
       .style('stroke-width', function(d) {
         return d._children ? 3 : 1;
       });
      }
      // Add labels for the nodes
      nodeEnter.append('text')
      .attr('dy', '.35em')
      .attr('x', function(d) {
        // Scale padding for label to the size of node
        var padding = (d.data.SizeOfNode || 10) + 3
        // PO
        // var padding = 10 + 3
        return d.children || d._children ? -1 * padding : padding;
      })
      .attr('text-anchor', function(d) {
        return d.children || d._children ? 'end' : 'start';
      })
      .style('font-size', options.fontSize + 'px')
      //.text(function(d) { return d.data.name; }); 
      // PO
      //.text(function(d) { return " "; });
      .text(function(d) { 

      		if(!options.pieNode){return d.data.name; } 
		isScar = JSON.parse(d.data.isScar);
		nodeLabel_sc = JSON.parse(options.nodeLabel_sc);
		if(isScar){return d.data.name;}
			if(nodeLabel_sc){ 
				if(d.data.parSize < options.nodeLabel_sc){ // TODO deprecate and adapt to show labels according to global layout 
		 			return d.data.ct;
			}
		}
		 return " ";
	 }); 

      // UPDATE
      var nodeUpdate = nodeEnter.merge(node);

      // Transition to the proper position for the node
      nodeUpdate.transition()
      .duration(duration)
      .attr('transform', function(d) {
        return 'translate(' + d.y + ',' + d.x + ')';
      });

      // Update the node attributes and style
      // PO choose depending on pie option
      // pie mode
      var cito = 1;
      if(options.pieNode ){
        nodeUpdate.select('g.node')
        .attr('outerRadius', function(d){
	    return p.size
	})
	// PO TODO verify action
	.style("fill", function(d) { 
         return d.color })
        .attr('cursor', 'pointer');
      }

      // PO circle mode
      if(!options.pieNode){
        nodeUpdate.select('circle.node')
        .attr('r', function(d) {
          return d.data.SizeOfNode || 10; // default radius is 10
        })
        .style('fill', function(d) {
          return d.data.fill || (d._children ? options.fill : '#fff');
        })
        .style('stroke-width', function(d) {
          return d._children ? 3 : 1;
        })
        .attr('cursor', 'pointer');
      }
      
      
      // Remove any exiting nodes
      var nodeExit = node.exit().transition()
      .duration(duration)
      .attr('transform', function(d) {
        return 'translate(' + source.y + ',' + source.x + ')';
      })
      .remove();

      // On exit reduce the node circles size to 0
      if(!options.pieNode){
      nodeExit.select('circle')
      .attr('r', 1e-6);
      }
      // On exit reduce the opacity of text labels
      nodeExit.select('text')
      .style('fill-opacity', 1e-6);

      // ****************** links section ***************************

      // Update the links...
      var link = svg.selectAll('path.link')
      .data(links, function(d) { return d.id; });

      // Enter any new links at the parent's previous position.
      var linkEnter = link.enter().insert('path', 'g')
      .attr('class', 'link')
      // Potentially, this may one day be mappable
      //.style('stroke-width', function(d) { return d.data.linkWidth || 1 })
      .style('stroke', function(d) {
        //if(options.pieNode)
      return (d.data.isScar ? "#555555" : "#999999")
      //return "#999999";
      })
      .style('stroke-width', function(d) {
        if(options.pieNode){
        	return (d.data.isScar ? 1.5 : 0.25)
	}
        return 1;
      })
      .style('stroke-dasharray', function(d) {
	if(d.data.Main != null){
	       	return (JSON.parse(d.data.Main) ? '1,0' : '2,3')
	}
      })
      .attr('d', function(d){
        var o = { x: source.x0, y: source.y0 }
        return diagonal(o, o)
      });

      // UPDATE
      var linkUpdate = linkEnter.merge(link);

      // Transition back to the parent element position
      linkUpdate.transition()
      .duration(duration)
      .attr('d', function(d){ return diagonal(d, d.parent) });

      // Remove any exiting links
      var linkExit = link.exit().transition()
      .duration(duration)
      .attr('d', function(d) {
        var o = {x: source.x, y: source.y}
        return diagonal(o, o)
      })
      .remove();

      // Store the old positions for transition.
      nodes.forEach(function(d){
        d.x0 = d.x;
        d.y0 = d.y;
      });

      // Creates a curved (diagonal) path from parent to the child nodes
      function diagonal(s, d) {

        path = 'M ' + s.y + ' ' + s.x + ' C ' +
        (s.y + d.y) / 2 + ' ' + s.x + ', ' +
        (s.y + d.y) / 2 + ' ' + d.x + ', ' +
        d.y + ' ' + d.x;

        return path
      }

      // Toggle children on click.
      function click(d) {
        if (JSON.parse(options.do_collapse)) {
          if (d.children) {
            d._children = d.children;
            d.children = null;
          } else {
            d.children = d._children;
            d._children = null;
          }
	}
        update(d);
        // Hide the tooltip after clicking
        tooltip.transition()
        .duration(100)
        .style('opacity', 0)
        // Update Shiny inputs, if applicable
        if (options.input) {
          var nest = {},
          obj = d;
          // Navigate up the list and recursively find parental nodes
          for (var n = d.depth; n > 0; n--) {
            nest[options.hierarchy[n-1]] = obj.data.name
            obj = obj.parent
          }
          Shiny.onInputChange(options.input, nest)
        }
      }

      // Show tooltip on mouseover
      function mouseover(d) {
        tooltip.transition()
        .duration(200)
        .style('opacity', .9);

        // Show either a constructed tooltip, or override with one from the data
        tooltip.html(
          d.data.tooltip || d.data.name + '<br>' +
          options.attribute + ': ' + d.data.WeightOfNode
        )
        // Make the tooltip font size just a little bit bigger
        .style('font-size', (options.fontSize + 1) + 'px')
        .style('left', (d3.event.layerX) + 'px')
        .style('top', (d3.event.layerY - 30) + 'px');
      }
      // Hide tooltip on mouseout
      function mouseout(d) {
        tooltip.transition()
        .duration(500)
        .style('opacity', 0);
      }
    }

    return {
      renderValue: function(x) {
        // Assigns parent, children, height, depth
        root = d3.hierarchy(x.data, function(d) { return d.children; });
        root.x0 = height / 2;
        root.y0 = 0;

        // Attach options as a property of the instance
        options = x.options;

        // Update the canvas with the new dimensions
        svg = svg.attr('transform', 'translate('
        + options.margin.left + ',' + options.margin.top + ')')

        // width and height, corrected for margins
        var heightMargin = height - options.margin.top - options.margin.bottom,
        widthMargin = width - options.margin.left - options.margin.right;
        // declares a tree layout and assigns the size
        treemap = d3.tree().size([heightMargin, widthMargin]);

        // Calculate a reasonable link length, if not otherwise specified
        if (!options.linkLength) {
          options.linkResponsive = true
          options.linkLength = widthMargin / options.hierarchy.length
          if (options.linkLength < 10) {
            options.linkLength = 10 // Offscreen or too short
          }
        }

        // Optionally collapse after the second level
        if (options.collapsed) root.children.forEach(collapse);
        update(root);

        // Collapse the node and all it's children
        function collapse(d) {
          if(d.children) {
            d._children = d.children
            d._children.forEach(collapse)
            d.children = null
          }
        }
      },

      resize: function(width, height) {
        // Resize the canvas
        d3.select(el).select('svg')
        .attr('width', width)
        .attr('height', height);

        // width and height, corrected for margins
        var heightMargin = height - options.margin.top - options.margin.bottom,
        widthMargin = width - options.margin.left - options.margin.right;

        // Calculate a reasonable link length, if not originally specified
        if (options.linkResponsive) {
          options.linkLength = widthMargin / options.hierarchy.length
          if (options.linkLength < 10) {
            options.linkLength = 10 // Offscreen or too short
          }
        }
        // Update the treemap to fit the new canvas size
        treemap = d3.tree().size([heightMargin, widthMargin]);
        update(root)
      },
      // Make the instance properties available as a property of the widget
      svg: svg,
      root: root,
      options: options
    };
  }
});

