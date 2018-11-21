// This piece of JS recieves a dict from the current pandas, and returns a csv download
var data = source;
var headers = Object.keys(data).join('\t').concat('\n');
for (var i = 0; i < data['Genes'].length; i++) {
    var rowarray = []
    for (var column in data){
        rowarray.push(data[column][i]);
    }
    var joined = rowarray.join("\t").concat('\n');
    headers = headers.concat(joined);
}
var filename = 'genes.csv';
var blob = new Blob([headers], { type: 'text/plain;charset=utf-8;' });
//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename);
} else {
    var link = document.createElement("a");
    link = document.createElement('a')
    link.href = URL.createObjectURL(blob);
    link.download = filename
    link.target = "_blank";
    link.style.visibility = 'hidden';
    link.dispatchEvent(new MouseEvent('click'))
}
