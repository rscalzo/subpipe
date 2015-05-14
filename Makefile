clean:
	find STAP Utils mydjango | grep '~$$' | xargs rm -rf
	find STAP Utils mydjango | grep '\.pyc$$' | xargs rm -rf

newref:
	rm -rf data/ref

newdata:
	rm -rf data/{sub,cand_info,mydjango}
	mkdir -p data/mydjango/{db,media}
	mkdir -p data/cand_info/{lcdata,skybot,xref}

newtestrun:
	rm -rf data/sub/2052 data/cand_info/*/3201
