// This piece of JS recieves a dict from the current pandas and a list of current Genes. Returns a fasta file download
var data = source;
var seqs = seqs;
var lines = "";
for (var i = 0; i < data['Genes'].length; i++) {
    header = ['>', data['Genes'][i]].join('').concat('\n');
    sequence = seqs[data['Genes'][i]].concat('\n');
    fasta = header.concat(sequence);
    lines = lines.concat(fasta);
}
var filename = 'genes.fasta';
var blob = new Blob([lines], { type: 'text/plain;charset=utf-8;' });
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
