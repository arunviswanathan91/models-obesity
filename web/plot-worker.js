// Python lives in a worker: tables and filters stay responsive during rendering.
let ready;
async function initialize(id){
  const runtimeURL=new URL('./runtime/',self.location.href).href;
  const {loadPyodide}=await import(runtimeURL+'pyodide.mjs');
  const py=await loadPyodide({indexURL:runtimeURL});
  postMessage({id,status:'Loading the plotting libraries…'});
  await py.loadPackage(['numpy','pandas','matplotlib','micropip']);
  py.globals.set("seaborn_wheel_url",runtimeURL+"seaborn-0.13.2-py3-none-any.whl");
  await py.runPythonAsync("import micropip\nawait micropip.install(seaborn_wheel_url, deps=False)");
  const response=await fetch('./plotting.py?v=manuscript-python-1');if(!response.ok)throw Error('Plotting code could not be loaded.');
  await py.runPythonAsync(await response.text());
  return py;
}
self.onmessage=async({data})=>{
  try{
    if(!ready){postMessage({id:data.id,status:'Starting Python for the first figure. This may take a minute…'});ready=initialize(data.id).catch(e=>{ready=null;throw e;});}
    const py=await ready;postMessage({id:data.id,status:'Drawing the selected results…'});
    py.globals.set('request_json',JSON.stringify(data.request));
    const result=await py.runPythonAsync('render(request_json)');
    postMessage({id:data.id,result:JSON.parse(result)});
  }catch(e){postMessage({id:data.id,error:String(e.message||e)});}
};
