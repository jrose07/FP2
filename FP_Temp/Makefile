MODE = None
NAME = vXXX
all:
ifneq ($(NAME), vXXX)
	@(find . -type f -name "*" -print0 | xargs -0 sed -i'' -e "s/vXXX/$(NAME)/g")
	@(mv vXXX/vXXX.tex vXXX/$(NAME).tex)
	@(mv vXXX $(NAME))
endif
	$(MAKE) -C $(NAME) MODE=$(MODE)
	cp $(NAME)/build/tex/$(NAME).pdf $(NAME)_rosenbaum_marschang.pdf

plots:
	$(MAKE) -C $(NAME) plot

clean:
	$(MAKE) -C $(NAME) clean

.PHONY: all clean
