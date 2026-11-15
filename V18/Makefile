MODE = None
NAME = v18
all:
ifneq ($(NAME), v18)
	@(find . -type f -name "*" -print0 | xargs -0 sed -i'' -e "s/v18/$(NAME)/g")
	@(mv v18/v18.tex v18/$(NAME).tex)
	@(mv v18 $(NAME))
endif
	$(MAKE) -C $(NAME) MODE=$(MODE)
	cp $(NAME)/build/tex/$(NAME).pdf $(NAME)_rosenbaum_marschang.pdf

plots:
	$(MAKE) -C $(NAME) plot

clean:
	$(MAKE) -C $(NAME) clean

.PHONY: all clean
