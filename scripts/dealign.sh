awk '{if (/^[^>]/) { gsub("-", ""); print } else { gsub(" ", ""); print }}' $1
