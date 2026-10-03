// Python lives in a worker: tables and filters stay responsive during rendering.
let ready;
async function initialize(){
  importScripts('https://cdn.jsdelivr.net/pyodide/v314.0.7/full/pyodide.js');
  const py=await loadPyodide({indexURL:'https://cdn.jsdelivr.net/pyodide/v314.0.7/full/'});
  postMessage({status:'Loading the plotting libraries…'});
  await py.loadPackage(['numpy','pandas','matplotlib','micropip']);
  await py.runPythonAsync("import micropip\nawait micropip.install('seaborn==0.13.2', deps=False)");
  const response=await fetch('./plotting.py');if(!response.ok)throw Error('Plotting code could not be loaded.');
  await py.runPythonAsync(await response.text());
  return py;
}
self.onmessage=async({data})=>{
  try{
    if(!ready){postMessage({status:'Starting Python for the first figure. This may take a minute…'});ready=initialize().catch(e=>{ready=null;throw e;});}
    const py=await ready;postMessage({status:'Drawing the selected results…'});
    py.globals.set('request_json',JSON.stringify(data.request));
    const result=await py.runPythonAsync('render(request_json)');
    postMessage({id:data.id,result:JSON.parse(result)});
  }catch(e){postMessage({id:data.id,error:String(e.message||e)});}
};
