

<!DOCTYPE html>

<html>

<head>
	<title>Leaflet Viz Demo</title>

	<meta charset="utf-8" />



	<meta name="viewport" content="width=device-width, initial-scale=1.0">



	<link rel="stylesheet" href="leaflet.css" />
     <link href="bootstrap.min.css" rel="stylesheet">
    
</head>


<jsp:useBean id="mt" scope="session" class="test.MotifShow" />
<jsp:setProperty name="mt" property="*" />


<body>

<center>
	<h1> Motif Visualization</h1>
	<br/>
	<div id="map" style="width: 500px; height: 500px"></div>
<br/><br/>
<table class="table" style="width: 500px;">
      <thead>
        <tr>
          <th>#</th>
          <th>Interval</th>
          <th>Length</th>
          <th>
			Show
		  </th>
        </tr>
      </thead>
      <tbody>
        <tr class="active">
          <th scope="row">1</th>
          <td>${mt.getRule(mt.get_start(),0)}</td>
          <td>${mt.getRule(mt.get_start()+0,1)}</td>
          <td>
		  	  		  		  <form name="frm1" method="post" action="MotifShow">
<input type="hidden" name="motif" />
<input id="b1" type="button" name="bt"  class="btn btn-info" style="width: 200px;" value="Show in Map" onclick="{document.frm1.motif.value=1;document.frm1.submit();}" />
</form>
</td>
 </tr>
			
        <tr>
          <th scope="row">2</th>
          <td>${mt.getRule(mt.get_start()+1,0)}</td>
          <td>${mt.getRule(mt.get_start()+1,1)}</td>
<td>
		  		  <form name="frm" method="post" action="MotifShow">
<input type="hidden" name="motif" />
<input id="b2" type="button" name="bt"  class="btn btn-info" style="width: 200px;" value="Show in Map" onclick="{document.frm.motif.value=2;document.frm.submit();}" />
</form>
</td>
        </tr>
        <tr class="active">
          <th scope="row">3</th>
          <td>${mt.getRule(mt.get_start()+2,0)}</td>
          <td>${mt.getRule(mt.get_start()+2,1)}</td>
		  <td>
		  		  		  <form name="frm3" method="post" action="MotifShow">
<input type="hidden" name="motif" />
<input id="b3" type="button" name="bt"  class="btn btn-info" style="width: 200px;" value="Show in Map" onclick="{document.frm3.motif.value=3;document.frm3.submit();}" />
</form>

		  </td>
        </tr>
        <tr>
          <th scope="row">4</th>
          <td>${mt.getRule(mt.get_start()+3,0)}</td>
          <td>${mt.getRule(mt.get_start()+3,1)}</td>
		  <td>
		  		  		  <form name="frm4" method="post" action="MotifShow">
<input type="hidden" name="motif" />
<input id="b4" type="button" name="bt"  class="btn btn-info" style="width: 200px;" value="Show in Map" onclick="{document.frm4.motif.value=4;document.frm4.submit();}" />
</form>
		  </td>
        </tr>
        <tr class="active">
          <th scope="row">5</th>
          <td>${mt.getRule(mt.get_start()+4,0)}</td>
          <td>${mt.getRule(mt.get_start()+4,1)}</td>
		  <td>
		  		  <form name="frm5" method="post" action="MotifShow">
<input type="hidden" name="motif" />
<input id="b5" type="button" name="bt"  class="btn btn-info" style="width: 200px;" value="Show in Map" onclick="{document.frm5.motif.value=5;document.frm5.submit();}" />
</form>

		  </td>
        </tr>
        <tr>
          <th scope="row">6</th>
          <td>${mt.getRule(mt.get_start()+5,0)}</td>
          <td>${mt.getRule(mt.get_start()+5,1)}</td>
		  <td>
		  		  <form name="frm6" method="post" action="MotifShow">
<input type="hidden" name="motif" />
<input id="b6" type="button" name="bt"  class="btn btn-info" style="width: 200px;" value="Show in Map" onclick="{document.frm6.motif.value=6;document.frm6.submit();}" />
</form>
		  </td>
        </tr>
      </tbody>
    </table>
	<table style="border: 10px">
	<tr>
	<td>
	<form name="n7r" method="post" action="MotifShow">
<input type="hidden" name="motif" />
<input id="b8" type="button" name="bt"  class="btn btn-primary" style="width: 200px;" value="Previous 6 Rules" onclick="{document.n7r.motif.value=8; document.n7r.submit();}" />
</form>
</td>
<td>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp </td>
<td>
<form name="n6r" method="post" action="MotifShow">
<input type="hidden" name="motif" />
<input id="b7" type="button" name="bt"  class="btn btn-primary" style="width: 200px;" value="Next 6 Rules" onclick="{document.n6r.motif.value=7;document.n6r.submit();}" />
</form>
</td>

<td>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp </td>
<td>
<form name="n8r" method="post" action="MotifShow">
<input type="hidden" name="motif" />
<input id="b7" type="button" name="bt"  class="btn btn-primary" style="width: 200px;" value="Return to Map View" onclick="{document.n8r.motif.value=9;document.n8r.submit();}" />
</form>
</td>


</tr>

</center>


	<script src="leaflet.js"></script>

	<script>

	
	
        
        var strMsg='${mt.getMsg(1)}';
        var arr = strMsg.split(' ');
        var i=0;
        
        for(i=0;i<arr.length; i++)
        	{
        		arr[i]=parseFloat(arr[i]);
        	}
        
        
        var newArr = [];
        while(arr.length) newArr.push(arr.splice(0,2));
        
        
		var map = L.map('map').setView([newArr[0][0],newArr[0][1]], 13);


		L.tileLayer('https://api.tiles.mapbox.com/v4/{id}/{z}/{x}/{y}.png?access_token=pk.eyJ1IjoibWFwYm94IiwiYSI6IjZjNmRjNzk3ZmE2MTcwOTEwMGY0MzU3YjUzOWFmNWZhIn0.Y8bhBaUMqFiPrDRW9hieoQ', {

			maxZoom: 18,

			attribution: 'Map data &copy; <a href="http://openstreetmap.org">OpenStreetMap</a> contributors, ' +

				'<a href="http://creativecommons.org/licenses/by-sa/2.0/">CC-BY-SA</a>, ' +

				'Imagery ?<a href="http://mapbox.com">Mapbox</a>',

			id: 'mapbox.streets'

		}).addTo(map);






var firstpolyline = new L.Polyline(newArr, {
color: 'blue',
weight: 3,
opacity: 0.5,
smoothFactor: 1

});
firstpolyline.addTo(map);


/*


L.circle(lat2[0], 1, {

	color: 'blue',

	fillColor: '#f03',

	fillOpacity: 0.5

}).addTo(map).bindPopup("I am step 1.");

L.circle(lat2[1], 1, {

	color: 'blue',

	fillColor: '#f03',

	fillOpacity: 0.5

}).addTo(map).bindPopup("I am step 2.");

L.circle(lat2[2], 1, {

	color: 'blue',

	fillColor: '#f03',

	fillOpacity: 0.5

}).addTo(map).bindPopup("I am step 3.");

		L.circle(lat2[3], 1, {

	color: 'blue',

	fillColor: '#f03',

	fillOpacity: 0.5

}).addTo(map).bindPopup("I am step 4.");


L.circle(lat2[4], 1, {

	color: 'blue',

	fillColor: '#f03',

	fillOpacity: 0.5

}).addTo(map).bindPopup("I am step 5.");

L.circle(lat2[5], 1, {

	color: 'blue',

	fillColor: '#f03',

	fillOpacity: 0.5

}).addTo(map).bindPopup("I am step 6.");


var point1 = new L.latLng(lat2[0][0], lat2[0][1]);
var point2 = new L.LatLng(lat2[1][0], lat2[1][1]);
var point3 = new L.LatLng(lat2[2][0], lat2[2][1]);
var point4 = new L.LatLng(lat2[3][0], lat2[3][1]);
var point5 = new L.LatLng(lat2[4][0], lat2[4][1]);
var point6 = new L.LatLng(lat2[5][0], lat2[5][1]);

var pointList = [point1, point2, point3, point4,point5,point6];




var firstpolyline = new L.Polyline(pointList, {
color: 'red',
weight: 3,
opacity: 0.5,
smoothFactor: 1

});
firstpolyline.addTo(map);


*/







//pointA = new L.LatLng(lat[1][0],lat[1][1]);
//pointB = new L.LatLng(lat[2][0],lat[2][1]);
//pointList = [pointA, pointB];

//firstpolyline = new L.polyline(lat, {
//color: 'red',
//weight: 3,
//opacity: 0.5,
//smoothFactor: 1

//});
//firstpolyline.addTo(map);


// zoom the map to the polyline
//map.fitBounds(polyline.getBounds());

	</script>

</body>

</html>
