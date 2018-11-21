// This piece of JS recieves a dict from the current pandas, and returns a csv download
//var data = source;
var a1 = ['Alex','Alexa','Tish'];
var a2 = [1000, 2000, 3000];
var a3 = [5,6,7];
var data = {'Genes': a1, 'income': a2, 'years_experience': a3};
var headers = Object.keys(data).join('\t').concat('\n');
print(headers);
for (var i = 0; i < data['Genes'].length; i++) {
    var rowarray = []
    for (var column in data){
        //print (column);
        rowarray.push(data[column][i]);
    }
    var joined = rowarray.join().concat('\n');
    headers = headers.concat(joined);
}
//for (const [ key, value ] of Object.entries(data)) {
   // do something with `key` and `value`
//}
//print(headers);
var filename = 'genes.csv';

