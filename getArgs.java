package inversion;
import java.util.Hashtable;

public class getArgs {
/*
	public static void main(String[] args) throws IOException, InterruptedException{
		parse(args);
		System.out.println(getString("input",false,""));
		System.out.println(getInt("int",false,100));
		System.out.println(getBoolean("boolean",false,false));
		System.out.println(getDouble("double",false,1.0));
	}
*/
	public static Hashtable<String, String> StringHash= new Hashtable<String, String>();
	public static String cmd="";
	public getArgs(String[] args){
		//parse args[i] to args[i+1] end at legnth-1
		for (int i=0;i<args.length-1;i++){
			if(args[i].length()>2 && args[i].substring(0,2).equals("--")){
				StringHash.put(args[i].substring(2), args[i+1]);
			}
		}
	}
	public void addUsage(String add){
		cmd=cmd.concat(add);
	}
	//input key,required parameter or not 1=required ,0=notRequired, default value
	public String getString(String key,boolean required,String a){
		if (StringHash.containsKey(key)){
			return StringHash.get(key);
		}else{
			if (required){
				System.out.println("No Key Name: "+key+"\n"+cmd);
				System.exit(1);
				return "";
			}else{
				return a;
			}
		}
	}
	public  boolean getBoolean(String key,boolean required,boolean a){
		if (StringHash.containsKey(key)){
			String tmp=StringHash.get(key);
			if (tmp.equalsIgnoreCase("true") || tmp.equalsIgnoreCase("t")){
				return true;
			}else if(tmp.equalsIgnoreCase("false")||tmp.equalsIgnoreCase("f")){
				return false;
			}else{
				System.out.println("Can't recognise boolean type for key: "+key+" . Please try True/T for true, False/F for false\n"+cmd);
				System.exit(1);
				return true;//useless
			}
		}else{
			if (required){
				System.out.println("No Key Name: "+key+"\n"+cmd);
				System.exit(1);
				return true;
			}else{
				return a;
			}
		}
	}
	public  int getInt(String key,boolean required,int a){
		if(StringHash.containsKey(key)){
			return Integer.parseInt( StringHash.get(key));
		}else{
			if (required){
				System.out.println("No Key Name: "+key+"\n"+cmd);
				System.exit(1);
				return 1;
			}else{
				return a;
			}
		}
	}
	public  double getDouble(String key,boolean required,double a){
		if(StringHash.containsKey(key)){
			return Double.parseDouble( StringHash.get(key));
		}else{
			if (required){
				System.out.println("No Key Name: "+key+"\n"+cmd);
				System.exit(1);
				return 1.0;
			}else{
				return a;
			}
		}	
	}
}
