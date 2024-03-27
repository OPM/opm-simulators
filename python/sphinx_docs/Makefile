.PHONY: docs

# Build the documentation locally for the current branch
# NOTE: You need to commit your changes before running this command
docs:
	@CURRENT_BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
	sphinx-versioned -m $$CURRENT_BRANCH -b $$CURRENT_BRANCH --git-root ../../