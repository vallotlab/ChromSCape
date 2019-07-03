# Rd
# description >> This function calculates the average values in a matrix per groups of lines or columns
# argument
# item >> mat >> Matrix to be averaged
# item >> margin >> Margin (1=lines, 2=columns)
# item >> groups >> List of groups
# item >> method >> mean or median
# value >> Averaged matrix
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.groupMat <- function(mat=NA,
						  margin=1,
						  groups=NA,
						  method="mean"
						  )
{
	if(!method %in% c("mean","median")){print("Method must be mean or median");break}
	if(!margin %in% 1:2){print("Margin must be 1 or 2");break}
	for(i in 1:length(groups))
	{
		if(margin==1){
			if(length(groups[[i]])==1){
				v <- mat[,groups[[i]]]
			}else{
				if(method=="mean"){v <- apply(mat[,groups[[i]]],margin,mean)}else{v <- apply(mat[,groups[[i]]],margin,median)}		
			}
			if(i==1){res <- matrix(v,ncol=1)}else{res <- cbind(res,v)}
		}else{
			if(length(groups[[i]])==1){
				v <- mat[groups[[i]],]
			}else{
				if(method=="mean"){v <- apply(mat[groups[[i]],],margin,mean)}else{v <- apply(mat[groups[[i]],],margin,median)}		
			}
			if(i==1){res <- matrix(v,nrow=1)}else{res <- rbind(res,v)}
		}
	}
	if(margin==1){rownames(res) <- rownames(mat);colnames(res) <- names(groups)}else{
		rownames(res) <- names(groups);colnames(res) <- colnames(mat)
	}
	res
}						  



