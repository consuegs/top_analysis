import getpass
import os
import configparser
import re


class BetterConfigParser(configparser.RawConfigParser):

    def get(self, section, option):
        result = configparser.RawConfigParser.get(self, section, option)
        result = self.__replaceSectionwideTemplates(result)
        return result

    def optionxform(self, optionstr):
        '''
        enable case sensitive options in .ini files
        '''
        return optionstr

    def __replaceSectionwideTemplates(self, data):
        '''
        replace <section|option> with get(section,option) recursively
        '''
        result = data
        findExpression = re.compile("((.*)\<(.*)\|(.*)\>(.*))*")
        groups = findExpression.search(data).groups()
        if not groups == (None, None, None, None, None): # expression not matched
            result = self.__replaceSectionwideTemplates(groups[1])
            result += self.get(groups[2], groups[3])
            result += self.__replaceSectionwideTemplates(groups[4])
        return result



def getPath(pathName, user=None):
    
    if user == None:
        user = getpass.getuser()
    
    if not os.path.isfile(os.path.join(os.path.dirname(__file__),"%s.ini"%user)):
        print("No path config found for user %s"%user)
        exit(1)
    
    config = BetterConfigParser()
    config.read(os.path.join(os.path.dirname(__file__),"%s.ini"%user))

    return config.get("paths",pathName)


if __name__ == "__main__":
    print(getPath("testPath"))
