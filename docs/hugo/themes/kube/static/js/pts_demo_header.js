
window.demoDescription = "In a field of points that revolves around a center, generate a voronoi cell diagram";

(function() {

  Pts.namespace( this );
  var space = new CanvasSpace("#pt").setup({bgcolor: "rgba(255,255,255,1)", resize: true, retina: true});
  var form = space.getForm();

  var pts = new Group();
  var timeOutId = -1;
  var header = null;


  space.add({ 

    // creatr 200 random points
    start:( bound ) => {
      pts = Create.distributeRandom( space.innerBound, 200 ); // Can we expand past innerBound?
      header = document.getElementById("header");
      // USE repel function from Delaunay example
    }, 

    animate: (time, ftime) => {
      // make a line and turn it into an "op" (see the guide on Op for more)
      let perpend = new Group( space.center.$subtract(0.1), space.pointer ).op( Line.perpendicularFromPt );
      pts.rotate2D( 0.0005, space.center );

      de = Create.delaunay(pts);
      triangles = de.delaunay();
      cells = de.voronoi();
      form.strokeOnly("#E8E8E8").polygons( cells );

      pts.forEach( (p, i) => {
        // for each point, find the perpendicular to the line
        let lp = perpend( p );
        var ratio = Math.min( 1, 1 - lp.$subtract(p).magnitude()/(space.size.x/2) );
        form.stroke(`rgba(255,255,255,${ratio}`, ratio*2).line( [ p, lp ] );
        form.fillOnly( ["#f03", "#09f", "#0c6"][i%3] ).point( p, 1.5, "circle" );
      });



      // header position
      if (header) {
        let top = window.pageYOffset || document.documentElement.scrollTop;
        let dp = top - space.size.y + 150;
        if (dp > 0) {
          header.style.top = `${dp * -1}px`;
        } else {
          header.style.top = "0px";
        }
      }

    },

    resize: () => {
      clearTimeout( timeOutId );
      setTimeout( () => {
        pts = Create.distributeRandom( space.innerBound, 200 );
      }, 500 );
    }

  });
  
  //// ----
  

  space.bindMouse().bindTouch().play();

})();;