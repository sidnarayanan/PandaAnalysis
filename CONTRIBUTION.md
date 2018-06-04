PandaAnalysis is a semi-collaborative project.
By this, I mean it was originally written by me, for me, to do my physics analyses.
However, it has proven useful to other people, and I welcome collaboration!
If you find a core bug in the code, please open a github issue to the relevant project, and I will look at it.
If you want to implement something for your use case, you can open a github issue, and depending on the complexity, I may be able to help you.
However, the general model for extending this package is pull request/review/merge. 

Here are the general rules for contribution:
- Your changes should not affect other analyses, unless explicitly agreed upon. This not only refers to the output information, but also the performance and memory usage.
- Test your code. I've resisted implementing CI because I have faith in you. Make sure it compiles and runs a reasonable test suite, which can be found in `PandaAnalysis/Flat/test/`
- Stick to the coding conventions of the project (spaces always, 2/4 for C/python, `if (` over `if(`, etc). Make sure everything is indented to the same level where appropriate.
- Code should be readable: no super long lines, use whitespace when it helps, etc.
- Absolutely no articles (definite or indefinite) in variable names :D. 
