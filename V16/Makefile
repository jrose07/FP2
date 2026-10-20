MODE = None
NAME = v16
all:
ifneq ($(NAME), v16)
	@(find . -type f -name "*" -print0 | xargs -0 sed -i'' -e "s/v16/$(NAME)/g")
	@(mv v16/v16.tex v16/$(NAME).tex)
	@(mv v16 $(NAME))
endif
	$(MAKE) -C $(NAME) MODE=$(MODE)
	cp $(NAME)/build/tex/$(NAME).pdf $(NAME)_rosenbaum_marschang.pdf

plots:
	$(MAKE) -C $(NAME) plot

clean:
	$(MAKE) -C $(NAME) clean

.PHONY: all clean
