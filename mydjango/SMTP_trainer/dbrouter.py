class dbrouter(object):
    """A router to control all database operations on models in
    the application"""

    def db_for_read(self, model, **hints):
        return self.__app_router(model)
    
    def db_for_write(self, model, **hints):
        return self.__app_router(model)
    
    def allow_relation(self, obj1, obj2, **hints):
        return obj1._meta.app_label == obj2._meta.app_label

    def allow_syncdb(self, db, model):
        return self.__app_router(model) == db
        
    def __app_router(self,model):
        if model._meta.app_label == 'candidate_new':
            return 'trainer_new'
        else:
            return 'default'
