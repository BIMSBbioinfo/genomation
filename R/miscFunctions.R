# you give a name of an RData object and a name of the variable and it assigns the given object to the variable
 Assigner = function(`_path`, `_name`){
        if(! is.character(`_path`) | !is.character(`_name`))
        stop('Both arguments should be characters!')
        load(`_path`)
        assign(`_name`, get(ls()[1]), parent.frame())
    }

# takes a name of the genome and loads it into a variable
# genome = GenomeLoader('BSgenome.Hsapiens.UCSC.hg19')
GenomeLoader = function(genome){
		
		require(genome, character.only=T)
		genome.name = unlist(strsplit(genome, split='\\.')) 
		return(get(genome.name[2]))		
	}
	
# given a vector and length smooths the vector to a given size
ScalerLarge = function(a, len, round.means=FALSE){
  
    if(length(a) < len)
		stop('vector can not be extended')
    s = unique(seq.int(1, length(a), length.out=len+1))
    starts = ceiling(s)[-length(s)]
    ends = floor(s)[-1]
    v = viewMeans(Views(a, starts, ends))
    if(round.means == TRUE){
      v=round(v)
    }
          
    return(v)
  }

TestScaler = function(){
  
  require(IRanges)
  error = character()
    
  # test1 
  r1 = Rle(c(1,2,3,4,5,5,5,3,2))
  l1 = 2
  if(!all(ScalerLarge(r1, l1) == c(3,4))){
    error = c(error, paste('Error:', length(r1), l1, "\n"))
  }
    
  # test2
  r2 = Rle(c(1,2,3,5,5,5,4,4,4))
  l2 = 3
  if(!all(ScalerLarge(r2, l2) == c(2,5,4))){
    error = c(error, paste('Error:', length(r2), l2, "\n"))
  }

  r3 = Rle(c(2,2,2,3,4,5,6,7,8))
  l3 = 4
  if(!all(ScalerLarge(r3, l3) == c(2,3,5,7))){
    error = c(error, paste('Error:', length(r3), l3, "\n"))
  }
 
  r4 = Rle(c(2,2,3,3,4,5,5,6,6))
  l4 = 5
  if(!all(ScalerLarge(r4, l4) == c(2,3,4,5,6))){
    error = c(error, paste('Error:', length(r4), l4, "\n"))
  }
 
  if(length(error) == 0){
    return(TRUE)
  }else{
    stop(error)
  } 
}
