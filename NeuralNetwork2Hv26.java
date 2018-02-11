import java.io.BufferedReader;         
import java.io.FileWriter; 
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Random;

// This code is modified from Microsoft MSDN BackPropagationProgram
// The default activation function is sigmoid
// Five steps to tanh activation function
// step 1 random [-0.5, 0.5]
// step 2 fun -> tanh
// step 3 label ->[-1 or 1]
// step 4 derive
// step 5 output
// tanh: 0.001	learning rate: ratel -- initial learning rate
// tanh: 0.4	momentum: rmoment -- initial momentum
// tanh: 0.2	b: alpha parameter in training function
// tanh: 2.0	dfac parameter for distance decay matrix

public class NeuralNetwork2Hv26{	
    int Inode; 
    int H1node; 
    int H2node;  
    int Onode;		
	
	float[][] wIH1;
	float[][] deltawIH1;	
	float[][] prevdwIH1;	
	float[][] wH1H2;
	float[][] deltawH1H2;
	float[][] prevdeltawH1H2;		
	float[][] wH2O;
	float[][] deltawH2O;
	float[][] prevdeltawH2O;	
	float[] inFeatures;		
	float[] inH1;
	float[] outH1;		
    float[] inH2;	
    float[] outH2;     	    
	float[] inO;			
	public float[] Ooutput;		
	float[] labelScore;		
	float learningRate;
	float momentum;	
	
	// tanh function
	// float b = 0.2f; This parameter is very important for tanh activiation function	
	public NeuralNetwork2Hv26(int ndim,int H1,int H2,int OnodeNum,float learningRate,float momentum){
	     this.Inode = ndim; 
	     this.H1node= H1; 
	     this.H2node = H2;  
	     this.Onode = OnodeNum;			     
	     this.wIH1 = new float[Inode][H1node];
	     this.prevdwIH1 = new float[Inode][H1node];		
	     this.wH1H2 = new float[H1node][H2node];
	     this.prevdeltawH1H2 = new float[H1node][H2node];		
	     this.wH2O = new float[H2node][Onode];
	     this.prevdeltawH2O = new float[H2node][Onode];		
	     this.inFeatures = new float[Inode];		
	     this.inH1 = new float[H1node];
	     this.outH1 = new float[H1node];		
	     this.inH2 = new float[H2node];	
	     this.outH2 = new float[H2node];     	    
	     this.inO = new float[Onode];			
	     this.Ooutput = new float[Onode];		
	     this.labelScore = new float[Onode];		
	     this.deltawIH1 = new float[Inode][H1node];		
	     this.deltawH1H2 = new float[H1node][H2node];	
	     this.deltawH2O = new float[H2node][Onode];
	     this.learningRate = learningRate;
	     this.momentum = momentum;	
	}
	
	public void loadWgt(String wgt) throws Exception{
		  resource t = new resource();
		  InputStream is=t.getClass().getResourceAsStream("/resource/" + wgt);   
		  BufferedReader br=new BufferedReader(new InputStreamReader(is));
		  int i,j;
		  for(i=0;i<Inode;i++){
		    for(j=0;j<H1node;j++){
		    	wIH1[i][j] = Double.valueOf(br.readLine()).floatValue();
			}
		  }
		  for(i=0;i<H1node;i++){
			for(j=0;j<H2node;j++){
				wH1H2[i][j] = Double.valueOf(br.readLine()).floatValue();
			}
		   }			
		   for(i=0;i<H2node;i++){
			 for(j=0;j<Onode;j++){
			    wH2O[i][j] = Double.valueOf(br.readLine()).floatValue();
			 }
		  }		   
		  br.close();
		  br=null;
	}
	
	public void saveWeights(String wgt) throws Exception{	    
		 int i, j;
		 FileWriter fw = new FileWriter(wgt);		 
	     for(i=0;i<Inode;i++){
			for(j=0;j<H1node;j++){
			   fw.append(wIH1[i][j] + "\n");
			}
	     }
		 for(i=0;i<H1node;i++){
		    for(j=0;j<H2node;j++){
			   fw.append(wH1H2[i][j]+"\n");
			}
		  }			
		  for(i=0;i<H2node;i++){
			for(j=0;j<Onode;j++){
			   fw.append(wH2O[i][j] + "\n");
			}
		  }		   
		  fw.flush();
		  fw.close();
		  fw = null;
	}

	public void updateWeights(){
		  int i,j;	  
		  float thisWeight,current,previous;
		  float thisDelta;
		  for(j=0;j<H1node;j++){
		    for(i=0;i<Inode;i++){  	  
		    	  thisWeight=wIH1[i][j]; 
		    	  current=deltawIH1[i][j]; 
		    	  previous=prevdwIH1[i][j];
		    	  thisDelta=this.learningRate*current + this.momentum*previous;
		    	  thisWeight+=thisDelta;	    	  
		    	  wIH1[i][j]=thisWeight;  
		    	  prevdwIH1[i][j]=thisDelta;	 	 
		    }
		  }		  
		  for(j=0;j<H2node;j++){
			   for(i=0;i<H1node;i++){  	  
			       thisWeight = wH1H2[i][j]; 
			       current = deltawH1H2[i][j]; 
			       previous = prevdeltawH1H2[i][j];
			       thisDelta = this.momentum*previous + this.learningRate*current;	    	  
			       thisWeight+=thisDelta;	    	  
			       wH1H2[i][j]=thisWeight;  
			       prevdeltawH1H2[i][j]=thisDelta;	 	 
			   }
		  }		  
		  for(j=0;j<Onode;j++){
		     for(i=0;i<H2node;i++){ 	  
		    	  thisWeight=wH2O[i][j]; 
		    	  current=deltawH2O[i][j]; 
		    	  previous=prevdeltawH2O[i][j];
		    	  thisDelta=this.momentum*previous + this.learningRate*current;
		    	  thisWeight+=thisDelta;	    	  
		    	  wH2O[i][j]=thisWeight;  
		    	  prevdeltawH2O[i][j]=thisDelta;	   	
		     }
		 }	
		 for(j=0;j<H1node;j++){
		    for(i=0;i<Inode;i++){
		    	deltawIH1[i][j]=0.0f;
		    }
		 }		 
		 for(j=0;j<H2node;j++){
			for(i=0;i<H1node;i++){
			    deltawH1H2[i][j]=0.0f;
			}
		 }	
		 for(j=0;j<Onode;j++){
		    for(i=0;i<H2node;i++){
		    	deltawH2O[i][j]=0.0f;
		    }
		 }	
	}
	
	public void randomizeWeights(){		  	  
		  int i,j;
		  Random r = new Random();
		  for(i=0;i<Inode;i++){
		    for(j=0;j<H1node;j++){
		    	wIH1[i][j]=(0.5f-r.nextFloat())*0.2f;
		    	deltawIH1[i][j]=0.0f;      
		    	prevdwIH1[i][j]=0.0f;	       
		    }
		  }		  
		  for(i=0;i<H1node;i++){
			 for(j=0;j<H2node;j++){
			    wH1H2[i][j]=(0.5f-r.nextFloat())*0.2f;
			    deltawH1H2[i][j]=0.0f;      
			    prevdeltawH1H2[i][j]=0.0f;	       
			  }
		  }	
		  for(i=0;i<H2node;i++){
			  for(j=0;j<Onode;j++){
			    wH2O[i][j]=(0.5f-r.nextFloat())*0.2f;	    	
			    deltawH2O[i][j]=0.0f;	    	      
			    prevdeltawH2O[i][j]=0.0f;	    	
			  }
		  }	
	} 
	
	public void forwardAndBackpropagation(){
		  this.forward(); 
		  this.backDerivation();
	}
	
	public void backDerivation(){
		  int i,j;		  
		  float[] dwO = new float[Onode];
		  float[] dwH2 = new float[H2node];		  
		  float[] dwH1 = new float[H1node];		    
		  for(i=0;i<Onode;i++)
			  dwO[i]=0;
		  for(i=0;i<H2node;i++) 
			  dwH2[i]=0.0f;		
		  for(i=0;i<H1node;i++) 
			  dwH1[i]=0.0f;		  
		  
		  // tanh: dwO[j]=b*((labelScore[j]-out[j])*(1.0f-out[j])*(1.0f+out[j]));
		  for(j=0;j<Onode;j++){
		    dwO[j]=((labelScore[j]-Ooutput[j])*Ooutput[j]*(1.0f-Ooutput[j]));  	   
		    for(i=0;i<H2node;i++){
		      deltawH2O[i][j]+=dwO[j]*outH2[i];
		      dwH2[i]+=dwO[j]*wH2O[i][j];
		    }		    
		  }	  			  
		  for(j=0;j<H2node;j++){
			dwH2[j]= dwH2[j]*outH2[j]*(1.0f-outH2[j]);
		    for(i=0;i<H1node;i++){
		      // dw1[i][j]+=dwH[j]*in2[i];
		      // dw1[i][j]+=dwH[j]*in2[i];
		      // tanh: dwH[j]= (float)(dwH[j]*b*(1.0f-out2[j])*(1.0f+out2[j]));
		      deltawH1H2[i][j]+=dwH2[j]*outH1[i];
		      dwH1[i]+=dwH2[j]*wH1H2[i][j];
		    }
		  }		  
		  for(j=0;j<H1node;j++){
			dwH1[j]= dwH1[j]*outH1[j]*(1.0f-outH1[j]);
			for(i=0;i<Inode;i++){
				deltawIH1[i][j]+=dwH1[j]*inFeatures[i];
			}
		  }	
	}
	
	public void forward(){
		  int i,j;		  
		  for(i=0;i<H1node;i++){
			inH1[i]=0.0f;
			for(j=0;j<Inode;j++){
			   inH1[i]+=wIH1[j][i]*inFeatures[j];
			}
			outH1[i]=1.0f/(1.0f+(float)Math.exp(-1.0f*inH1[i]));
		  }		  
		  for(i=0;i<H2node;i++){
			inH2[i]=0.0f;
		    for(j=0;j<H1node;j++){
		       inH2[i]+=wH1H2[j][i]*outH1[j];
		    }
		    outH2[i]=1.0f/(1.0f+(float)Math.exp(-1.0f*inH2[i]));
		  }
		  
		  for(i=0;i<Onode;i++){
			inO[i]=0.0f;
		    for(j=0;j<H2node;j++){
		       inO[i]+=wH2O[j][i]*outH2[j];
		    }
		    Ooutput[i]=1.0f/(1.0f+(float)Math.exp(-1.0*inO[i]));
		  }
	}		
	
	public void solver(float[] score,float[] label){
		  this.feedFeatures(score, label);
		  this.forwardAndBackpropagation();
		  this.updateWeights();	   
	}
	
	public void feedFeatures(float[] score,float[] label){
		  int i;		 
		  for(i=0;i<this.Inode;i++){
			  this.inFeatures[i] = score[i];
		  }		  
		  for(i=0;i<label.length;i++){
			  this.labelScore[i] = label[i];
		  }
	}	
	
	// This is the activiation function f(x) = 1/(1+exp(alpha*x))
	// The derivative f'(x) = alpha*f(x)*(1-f(x))			
	// private float tfunc(float xoutp){			   		
	//	return (float)Math.tanh(alpha * xoutp);
	// }
	
}