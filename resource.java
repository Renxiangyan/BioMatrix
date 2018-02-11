
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;


public class resource {

	public String getResource(String myfile) throws IOException{     
		         URL fileURL=this.getClass().getResource(myfile);   
		         return fileURL.getFile();  
    }  
	
	public static void main(String[] args) throws Exception, IOException {			
		
		resource t = new resource();
		InputStream is=t.getClass().getResourceAsStream("/resource/r");   
		BufferedReader br=new BufferedReader(new InputStreamReader(is));  
		String line = "";
		while((line = br.readLine()) != null){	
			line = line.replaceFirst("\\s+", "");
			System.out.println(line);
		}
		br.close();
	}
	
}






















